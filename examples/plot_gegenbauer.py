# %%
from scipy.special import gegenbauer
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
figname = 'First_few_Gegenbauer_polynomials'
figtitle = r'First few Gegenbauer polynomials $C_{\lambda}^{\nu}(t)$ for $\lambda \in [0 \, ..\, 3]$, $l \in [1 \, ..\, 4]$'
savedir = '../data/' + figname
# %%
z_grid = np.linspace(-1, 1, 10000)
# create array of axes in one figure
fig, ax = plt.subplots(2, 2, figsize=(12, 12), constrained_layout=True)
fig.suptitle(figtitle, fontsize=22)
row_index = 0
col_index = 0

# loop over n and l
for lamb in range(0, 4):
    plt.sca(ax[row_index, col_index])

    for nu in range(1, 5):
        C_vals = gegenbauer(lamb, nu)(z_grid)
    
        plt.plot(z_grid, C_vals, label=r'$\lambda = {0}§´$, $\nu = {1}$'.format(lamb, nu))

    plt.legend()
    plt.ylim(-10, 10)

    # manage labeling only the left-most y axes and the bottom-most x axes
    if col_index == 0:
        plt.ylabel(r'$C_{\lambda}^\nu(t)$')

    if row_index == 1:
        plt.xlabel('$t$')


    # manage placing the axes at the correct position
    row_index += 1
    
    if row_index >= 2:
        col_index += 1
        row_index = 0

plt.savefig(fname=savedir, dpi = 300)
plt.show()

# %%

# %%
