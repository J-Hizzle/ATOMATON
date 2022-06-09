# %%
from scipy.special import sph_harm
from scipy.constants import pi
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
figname = 'First_few_spherical_harmonics'
figtitle = 'First few spherical harmonics $Y_{l, m}(\theta, \phi)$ for $\theta \in [0, \pi]$, $\phi \in [0, 2 \pi]$'
savedir = '../data/' + figname
# %%
# create array of axes in one figure
fig, ax = plt.subplots(2, 2, figsize=(12, 12), constrained_layout=True)
fig.suptitle(figtitle, fontsize=22)
row_index = 0
col_index = 0

# loop over n and l
for n in range(1, 5):
    plt.sca(ax[row_index, col_index])

    # manage labeling only the left-most y axes and the bottom-most x axes
    if col_index == 0:
        plt.ylabel('$Y_{l, m}(\theta, \phi)$')
        #plt.xlim(0, 4)
    elif col_index == 1:
        #plt.xlim(0, 2)

    if row_index == 1:
        plt.xlabel('$y/n$')


    y_grid = p_grid * n/Z

    for l in range(0, n):
        G_vals = radial_function(n, l, Z, y_grid)
    
        plt.plot(y_grid/n, G_vals * np.sqrt(Z), label='$n$ = {0}, $l$ = {1}'.format(n, l))

    plt.legend()
    plt.ylim(-3.5, 3.5)

    # manage placing the axes at the correct position
    row_index += 1
    
    if row_index >= 2:
        col_index += 1
        row_index = 0

plt.savefig(fname=savedir, dpi = 300)
plt.show()
