# %%
from atomaton.momentum_space import calc_momentum_distribution_function
from atomaton.momentum_space_rel_fromtaluk import calc_momentum_distribution_function_rel
from atomaton.auxiliary_functions import rescale_axes
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
Z_vals = [1, 80, 137]
n_vals = [1, 2, 3]
l_vals = [0, 0, 0]
kappa_vals = [-1, -1, -1]
p_grid = np.linspace(0, 10, 100)

figname = 'comparison_momentum_distribut4ion'
figtitle = 'Comparison between Schrödinger and Dirac momentum distribution functions'
savedir = '../data/' + figname
# %%
# create array of axes in one figure
fig, ax = plt.subplots(3, 3, figsize=(12, 12), constrained_layout=True, gridspec_kw={'width_ratios': [1, 1, 1]})
fig.suptitle(figtitle, fontsize=22)

for row_index in range(0, 3):
    for col_index in range(0, 3):
        # pick correct set of quantum numbers
        n_val = n_vals[col_index]
        l_val = l_vals[col_index]
        kappa_val = kappa_vals[col_index]
        Z_val = Z_vals[row_index]

        # pick the indexed axis as the current plot
        plt.sca(ax[row_index, col_index])

        # construct and evaluate Schrödinger momentum distribution function
        dist_nr = calc_momentum_distribution_function(n_val, l_val, Z_val, p_grid)

        # construct and evaluate Dirac momentum distribution function
        dist_rel = calc_momentum_distribution_function_rel(n_val, kappa_val, Z_val, p_grid, alpha_val=0)

        print(type(dist_rel))

        # rescale axes
        p_grid_rescaled, dist_nr_rescaled, rescale_index = rescale_axes(p_grid, dist_nr, isovalue=0.999)

        dist_rel_rescaled = dist_rel[:rescale_index]

        plt.plot(p_grid_rescaled, dist_nr_rescaled)
        plt.plot(p_grid_rescaled, dist_rel_rescaled)

        plt.yscale('log')

plt.show()
# %%
plt.show()
# %%
