import argparse
import numpy as np
import os
from subprocess import run, PIPE
import time


#ceu simulation:
#for neutral:
#  mutation rate ~U(1.25e-08, 2.5e-08)
#      mutrates = runif(1, 1.25e-8, 2.5e-8);
#  recombination rate ~U(1.25e-08, 2.5e-08)
#      recombrates = runif(1, 1.25e-8, 2.5e-8);
#
#
#for hard sweep
#  selection coef ~U(0.0001, 0.02)
#      selcoef = runif(1, 0.0001, 0.02);
#  segregating freq of site under selection ~U(0.01, 0.99)
#      segfreq = runif(1, 0.01, 0.99);


parser = argparse.ArgumentParser()
parser.add_argument('--outdir', default = 'data')
parser.add_argument('--array', required=True)
parser.add_argument('--jobs', required=True)
parser.add_argument('--jobsstart', default=101)
args = parser.parse_args()

outdir = args.outdir
array = args.array
jobs = int(args.jobs)
jobsstart = int(args.jobsstart)

SLIM = '/scratch1/bldinh/csci567/slimbuilds/build/slim'
RELATEFILEFORMATS = '/project/chia657_28/programs/relate_v1.2.0_x86_64_static/bin/RelateFileFormats'
RELATE = '/project/chia657_28/programs/relate_v1.2.0_x86_64_static/bin/Relate'


def write_rate_to_file(rate, fp):
    with open(fp, 'w') as g:
        g.write(f'{rate}\n')


def run_hardsweep_slim_command(rep, m_rate, r_rate, s_coef):
    write_rate_to_file(m_rate, 'sim.mrate')
    write_rate_to_file(r_rate, 'sim.rrate')
    write_rate_to_file(s_coef, 'sim.scoef')

    template = ("""// Sweep simulation with tree sequence recording
// need to define the following Eidos constants:
// `paramF`, `outPref`, `selcoef`, `mutgen`, `min_AF`, `max_AF`

initialize() {
    paramF = '/scratch1/bldinh/csci567/sim_ceu/dem_eg.param';

    defineConstant("outPref", "simhardsweep");
    defineConstant("mutgen", 500);
    defineConstant("min_AF", 0.2);
    defineConstant("max_AF", 0.95);
    params = readFile(paramF);

    gen = c();
    Ne = c();

    for (l in 2:(size(params)-1)){
        gen_Ne = asInteger(strsplit(params[l]));
        gen = c(gen, gen_Ne[0]);
        Ne = c(Ne, gen_Ne[1]);
    }

    defineConstant("GOI", gen[1:(size(gen)-2)]);
    defineConstant("Ne_GOI", Ne[1:(size(Ne)-2)]);

    defineConstant("N_0", Ne[0]); // N_0
    defineConstant("last", gen[size(gen)-1]);  // end gen

    initializeTreeSeq();
    """
    "initializeMutationRate(" + f"{m_rate}"+ ");\n"
    "initializeRecombinationRate(" + f"{r_rate}" + ");\n"

    """
    initializeMutationType("m1", 0.5, "f", 0.0);
    """

    "initializeMutationType(\"m2\", 0.5, \"f\"," + f"{s_coef}" + ");\n"
    """
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, 1e5-1);
}

1 {
    defineConstant("simID", getSeed());
    sim.addSubpop("p1", N_0);
    sim.rescheduleScriptBlock(s1, generations=GOI);
    sim.rescheduleScriptBlock(s2, mutgen-100, mutgen-100);
    sim.rescheduleScriptBlock(s3, mutgen, mutgen);
    sim.rescheduleScriptBlock(s4, last, last);
}

s1 2 {
    if (sum(GOI==sim.generation)){
        p1.setSubpopulationSize(Ne_GOI[GOI==sim.generation]);
    }
}

s2 3 late() {
    sim.treeSeqOutput("/scratch1/bldinh/csci567/sim_ceu/ceu_hardsweep/tmp/slim_" + simID + ".trees");
}

s3 3 late() {
    target = sample(p1.genomes, 1);
    target.addNewDrawnMutation(m2, 5e4);
}

1: late() {
    if (sim.generation > mutgen & sim.countOfMutationsOfType(m2) == 0){
        fixed = (sum(sim.substitutions.mutationType == m2) == 1);

        //if (fixed){
        //    cat(simID + ": FIXED - RESTARTING\\n");
        //}
        //else{
        //    cat(simID + ": LOST - RESTARTING\\n");
        //}
        // go back to generation `mutgen-100`
        sim.readFromPopulationFile("/scratch1/bldinh/csci567/sim_ceu/ceu_hardsweep/tmp/slim_" + simID + ".trees");

        // start a newly seeded run
        setSeed(rdunif(1, 0, asInteger(2^32) - 1));
    }
}

s4 4 late() {
    if (sim.countOfMutationsOfType(m2) == 0){
        mut_freq = -1; // either fixed or lost
    } else {
        mut_freq = sim.mutationFrequencies(NULL, sim.mutationsOfType(m2));
    }
    if (mut_freq < min_AF | mut_freq > max_AF){
        //cat(simID + ": FREQ OOB - RESTARTING\\n");

        // go back to generation `mutgen-100`
        sim.readFromPopulationFile("/scratch1/bldinh/csci567/sim_ceu/ceu_hardsweep/tmp/slim_" + simID + ".trees");

        // start a newly seeded run
        setSeed(rdunif(1, 0, asInteger(2^32) - 1));
        return;
    }

    sim.treeSeqOutput(paste(c(outPref, ".trees"), sep=""));
    //writeFile(paste(c(outPref,'.selcoef'), sep=""), asString(selcoef));
    //writeFile(paste(c(outPref,'.daf'), sep=""), asString(mut_freq));

    //p1.outputMSSample(198, replace=F, filePath=paste(c(outPref,".ms"), sep=""));
    p1.outputVCFSample(99, replace=F, outputMultiallelics=F, filePath=paste(c(outPref, ".trees.vcf"), sep=""));

    cat(c("%%", mutgen, mut_freq, outPref, "\\n"), sep='\\t'); // make sure the treeseq file is saved before printing meta-data
    deleteFile("/scratch1/bldinh/csci567/sim_ceu/ceu_hardsweep/tmp/slim_" + simID + ".trees");

    sim.simulationFinished();
}
""")

    cmd = f'{SLIM}'
    p = run(cmd.split(), stdout=PIPE, input=template, encoding='ascii')
    return p.returncode


starting_dir = os.getcwd()
print(starting_dir)

for replicate in range(jobsstart,jobsstart+jobs):

    #
    # Make the output directory
    #
    outprefix = f'{starting_dir}/{outdir}/cpu{array}/rep{replicate}'
    cmd = f'mkdir -p {outprefix}'
    run(cmd.split())
    os.chdir(outprefix)
    time.sleep(1)

    #
    # SLiM
    #
    mut_rate = np.random.uniform(low=1.25e-8, high=2.5e-8)
    rec_rate = np.random.uniform(low=1.25e-8, high=2.5e-8)
    sel_coef = np.random.uniform(low=0.0001, high=0.02)
    print(f'rep{replicate}: mut:{mut_rate}, rec:{rec_rate}, sel:{sel_coef}')

    code = run_hardsweep_slim_command(replicate, mut_rate, rec_rate, sel_coef)
    print(f'rep{replicate} returned {code}')


    #
    # Relate
    #
    cmd = f'{RELATEFILEFORMATS} --mode ConvertFromVcf --haps simhardsweep.trees.haps --sample simhardsweep.trees.sample -i simhardsweep.trees'
    #print(cmd)
    run(cmd.split())

    cmd = f'{RELATE} --mode All -m {mut_rate} -N 5002870 --haps simhardsweep.trees.haps --sample simhardsweep.trees.sample -o simrelate --map /scratch1/bldinh/csci567/sim_ceu/dummy.uniform.map'
    #print(cmd)
    run(cmd.split())
    print()



