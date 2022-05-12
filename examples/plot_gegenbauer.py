# %%
from scipy.special import gegenbauer
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
z_grid = np.linspace(-1, 1, 10000)
# %%
figname = 'First_few_Gegenbauer_polynomials'
figtitle = 'First few Gegenbauer polynomials $C_{n_r}^{l + 1}(z)$ for $z \in [-1, 1]$'
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

    for l in range(0, n):
        C_vals = gegenbauer(n - l - 1, l + 1)(z_grid)
    
        plt.plot(z_grid, C_vals, label='$n$ = {0}, $l$ = {1}'.format(n, l))

    plt.legend()
    plt.ylim(-10, 10)

    # manage labeling only the left-most y axes and the bottom-most x axes
    if col_index == 0:
        plt.ylabel('$C_{n_r}^{l + 1}(z)$')

    if row_index == 1:
        plt.xlabel('$z$')


    # manage placing the axes at the correct position
    row_index += 1
    
    if row_index >= 2:
        col_index += 1
        row_index = 0

plt.savefig(fname=savedir, dpi = 300)
plt.show()

# %%

# %%
