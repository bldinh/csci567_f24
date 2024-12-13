import numpy as np

fp = 'all.train.x.y.fixed.txt'

train_x = np.full((2000000,600), np.inf)
train_y = []
idx = 0
print('parsing file')
with open(fp) as f:
    for rline in f:
        line = [int(v) for v in rline.split()]
        train_x[idx] = np.array(line[:-1])
        train_y.append(line[-1])
        idx += 1

print(f'writing y file: {len(train_y)} lines')
with open('all.train.fixed.y.txt', 'w') as g:
    for y in train_y:
        g.write(f'{y}\n')

print('calculating mean and sigma')
mu = np.mean(train_x, axis=0)
sigma = np.std(train_x, axis=0)

print('calculating centered x')
train_x_centered = (train_x - mu)/sigma
print('writing centered x file')
train_x_centered.tofile('all.train.fixed.x.txt')

