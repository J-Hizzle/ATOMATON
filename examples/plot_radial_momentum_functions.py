# %%
from atomaton.momentum_space import calc_real_radial_momentum_function
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
Z = 1
p_grid = np.linspace(0, 10, 10000)

figname = 'First_few_radial_momentum_functions_nocap'
figtitle = 'First few radial functions in momentum space $\sqrt{Z} G_{n, l}(p)$ for $n \in [1 \, ..\, 4]$, $l \in [0 \, ..\, (n - 1)]$'
savedir = '../data/' + figname
# %%
# create array of axes in one figure
fig, ax = plt.subplots(2, 2, figsize=(12, 12), constrained_layout=True, gridspec_kw={'width_ratios': [2, 1]})
#fig.suptitle(figtitle, fontsize=22)
row_index = 0
col_index = 0

# loop over n and l
for n in range(1, 5):
    plt.sca(ax[row_index, col_index])

    # manage labeling only the left-most y axes and the bottom-most x axes
    if col_index == 0:
    #    plt.ylabel(r'$\sqrt{Z} G_{n, l}(p) \hspace{5mm} (a.u.)$')
        plt.xlim(0, 4)
    elif col_index == 1:
        plt.xlim(0, 2)

    #if row_index == 1:
    #    plt.xlabel('$p/Z \hspace{5mm} (a.u.)$')

    for l in range(0, n):
        G_vals = calc_real_radial_momentum_function(n, l, Z, p_grid)
        plt.plot(p_grid/Z, G_vals * np.sqrt(Z), label='$n$ = {0}, $l$ = {1}'.format(n, l))

    plt.legend()
    plt.ylim(-3.5, 3.5)

    # manage placing the axes at the correct position
    row_index += 1
    
    if row_index >= 2:
        col_index += 1
        row_index = 0

# label x and y axes for entire figure
fig.supylabel(r'$\sqrt{Z} G_{n, l}(p) \hspace{5mm} (a.u.)$', fontsize=25)
fig.supxlabel('$p/Z \hspace{5mm} (a.u.)$', fontsize=25)


plt.savefig(fname=savedir, dpi = 300)
plt.show()
# %%
