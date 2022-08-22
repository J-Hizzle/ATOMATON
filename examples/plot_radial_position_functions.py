# %%
from atomaton.position_space import calc_radial_position_function
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
Z = 1
r_grid = np.linspace(0, 100, 10000)

figname = 'First_few_radial_position_functions'
figtitle = r'First few radial functions $P_{n, l}(r)/\sqrt{Z}$ for $n \in [1 \, .. \, 4]$, $l \in [0 \, .. \, (n - 1)]$'
savedir = '../data/' + figname


# %%
# create array of axes in one figure
fig, ax = plt.subplots(2, 2, figsize=(12, 12), constrained_layout=True, gridspec_kw={'width_ratios': [1, 2]})
fig.suptitle(figtitle, fontsize=22)
row_index = 0
col_index = 0

# loop over n and l
for n in range(1, 5):
    plt.sca(ax[row_index, col_index])

    # manage labeling only the left-most y axes and the bottom-most x axes
    if col_index == 0:
        plt.ylabel(r'$P_{n, l}(r)/\sqrt{Z}\hspace{5mm} (a.u.)$')
        plt.xlim(0, 20)
    elif col_index == 1:
        plt.xlim(0, 40)

    if row_index == 1:
        plt.xlabel('$Zr \hspace{5mm} (a.u.)$')

    for l in range(0, n):
        P_vals = calc_radial_position_function(n, l, Z, r_grid)
    
        plt.plot(r_grid * Z, P_vals/np.sqrt(Z), label='$n$ = {0}, $l$ = {1}'.format(n, l))

    plt.legend()
    plt.ylim(-0.8, 0.8)

    # manage placing the axes at the correct position
    row_index += 1
    
    if row_index >= 2:
        col_index += 1
        row_index = 0

plt.savefig(fname=savedir, dpi = 300)
plt.show()
# %%
