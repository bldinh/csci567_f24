from collections import deque
import numpy as np
import argparse
import random
import sys
import tskit


parser = argparse.ArgumentParser()
parser.add_argument('--array', default=1)
parser.add_argument('--replicates', default=1)
parser.add_argument('--outdir', default='new_data')
parser.add_argument('--samplehap', default=198)

args = parser.parse_args()

array = int(args.array)
num_iterations = int(args.replicates)
haplotypes_to_sample = int(args.samplehap)
outdir = args.outdir

all_tree_file = 'current.500.paths.trees.txt'
with open(all_tree_file) as f:
    all_tree_fp = f.read().splitlines()

list_of_paths_to_trees = []
for idx, fp in enumerate(all_tree_fp):
    if idx % 100 == array:
        list_of_paths_to_trees.append(fp)

print(f'total tree paths: {len(all_tree_fp)}')
print(f'tree paths for current worker:: {len(list_of_paths_to_trees)}')

#list_of_paths_to_trees = list_of_paths_to_trees[:1]


def extract_feature_vector(list_of_trees, s_coef, list_of_time_points, treeseq):
    #take in 5 trees (2 left, tree of interest, 2 right)
    #return list of features and class based on selection coef

    #logic based on neutral or sweep:
    if s_coef > 0:
        sweep_class = '1'
        sweep = True
    else:
        sweep_class = '0'
        sweep = False

    f_vector = []
    f_vector_prefix = [list_of_trees[0].num_lineages(t) for t in list_of_time_points] + [list_of_trees[1].num_lineages(t) for t in list_of_time_points]
    f_vector_suffix = [list_of_trees[3].num_lineages(t) for t in list_of_time_points] + [list_of_trees[4].num_lineages(t) for t in list_of_time_points]

    site = ''
    if sweep:
        for s in list_of_trees[2].sites():
            if s.position == 50000:
                site = s
    else:
        #randomly choose
        #make list of sites that intersect with variant sites...
        valid_sites = []

        for s in list(list_of_trees[2].sites()):
            for v in treeseq.variants():
                if v.site.position == s.position and len(s.mutations) == 1:
                    valid_sites.append(s)
        #subset from that
        #repeating logic from neutral tree selection
        #site = random.sample(list(list_of_trees[2].sites()),1)[0]
        site = random.sample(valid_sites,1)[0]


    d_af = -1
    if len(site.mutations) > 1:
        print('need logic for multiple mutations at the position')
    elif len(site.mutations) == 1:
        #want to grab all nodes below this node
        mut_node_id = site.mutations[0].node
        node_ids = [mut_node_id]
        queue = deque(node_ids)
        while queue:
            node = queue.popleft()
            for child in list_of_trees[2].children(node):
                node_ids.append(child)
                queue.append(child)
        #print(f'num of children nodes under mut node: {len(node_ids)}') #e.g. 360+

        mut_ts = ts_samp.subset(node_ids) #grab the tree sequence that only has mut. nodes
        der_lin_counts = []
        anc_lin_counts = []
        for mut_tree in mut_ts.trees():
            mut_left, mut_right = mut_tree.interval
            if not mut_left < site.position or not mut_right > site.position:
                continue
            for t in list_of_time_points:
                total_lineages = list_of_trees[2].num_lineages(t)
                deriv_lineages = mut_tree.num_lineages(t)
                der_lin_counts.append(deriv_lineages)
                anc_lineages = total_lineages - deriv_lineages
                anc_lin_counts.append(anc_lineages)
                #print('\t'.join([str(t) for t in [t, total_lineages, deriv_lineages, anc_lineages]]))
        f_vector = f_vector_prefix + der_lin_counts + anc_lin_counts + f_vector_suffix

        for v in treeseq.variants():
            if v.site.position == site.position:
                #print(v.genotypes)
                #print(len(v.genotypes))
                an = sum([1 for gt in v.genotypes if gt == 1 or gt == 0])
                ac = sum([1 for gt in v.genotypes if gt == 1])
                d_af = ac/an
                print(f'{len(f_vector)}: daf - {ac}/{an}={d_af:.4f}')
                break


    return (f_vector, sweep_class, d_af)



out_x = f'{outdir}/job{array}.features.class.txt'
out_daf = f'{outdir}/job{array}.daf.scoef.txt'
out_errors = f'{outdir}/job{array}.errors.txt'
with open(out_x, 'w') as g:
    with open(out_daf, 'w') as h:
        with open(out_errors, 'w') as i:
            for fpidx, fp in enumerate(list_of_paths_to_trees):
                print(f'fpidx: {fpidx}')

                selcoef_fp = fp.replace('simhardsweep.trees','sim.scoef')
                with open(selcoef_fp) as f:
                    selcoef = float(f.read().splitlines()[0])

                ts = tskit.load(fp)
                #print(ts)

                all_samp = ts.samples()

                for replicate in range(1,num_iterations+1):
                    print(f'rep{replicate}')
                    samp = list(random.sample(list(all_samp), haplotypes_to_sample))
                    #print(f'num samp: {len(samp)}')

                    ts_samp = ts.simplify(samples=samp)
                    time_pts = np.logspace(0, np.log10(ts_samp.max_root_time-1), num=100)
                    trees = ts_samp.aslist()


                    #randomly sample an index for the neutral location
                    print('start while loop')
                    max_iter = 5000
                    cur_iter = 0
                    #TODO logic to break/resample go to next iter
                    while True:
                        middle_idx = random.sample(range(2,len(trees)-2),1)[0]
                        tree_idxs = list(range(middle_idx-2, middle_idx+3))
                        if all([(trees[i].interval[1] < 50000) or (trees[i].interval[0] > 50000) for i in tree_idxs]) and (len(list(trees[middle_idx].sites())) > 0):
                            #print(list(trees[middle_idx].sites()))
                            matching_pos = []
                            for s in trees[middle_idx].sites():
                                for v in ts_samp.variants():
                                    if v.site.position == s.position and len(s.mutations) == 1:
                                        matching_pos.append(s.position)
                            if len(matching_pos) == 0:
                                continue
                            else:
                                break
                        cur_iter += 1
                    print('out of while loop')

                    current_list_of_trees = [trees[middle_idx-2], trees[middle_idx-1], trees[middle_idx], trees[middle_idx+1], trees[middle_idx+2]]
                    feature_vector, y_class, daf = extract_feature_vector(current_list_of_trees, 0, time_pts, ts_samp)

                    if daf >= 0:
                        g.write(' '.join([str(v) for v in feature_vector])+f' {y_class}\n')
                        h.write(f'{daf}\t0.0\n')
                    else:
                        i.write(f'no variant site matching mutation site for tree {trees[middle_idx].interval[0]} - {trees[middle_idx].interval[1]} in {fp}')

                    # mutation location
                    for idx, tree in enumerate(trees[:-2]):
                        left, right = tree.interval
                        #randomly select neutral mutation, check that it doesn't overlap sweep location

                        #find the hardsweep location
                        if not left < 50000 or not right > 50000:
                            pass
                        else:
                            current_list_of_trees = [trees[idx-2], trees[idx-1], trees[idx], trees[idx+1], trees[idx+2]]
                            feature_vector, y_class, daf = extract_feature_vector(current_list_of_trees, selcoef, time_pts, ts_samp)

                            if daf >= 0:
                                g.write(' '.join([str(v) for v in feature_vector])+f' {y_class}\n')
                                h.write(f'{daf}\t{selcoef}\n')
                            else:
                                i.write(f'no variant site matching mutation site for tree {left} - {right} in {fp}')
    print()























