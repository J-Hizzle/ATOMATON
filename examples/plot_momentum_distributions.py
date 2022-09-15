# %%
from atomaton.momentum_space import calc_momentum_distribution_function
from atomaton.momentum_space_rel_num import calc_momentum_distribution_function_rel
from atomaton.auxiliary_functions import rescale_axes
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
Z_vals = [1, 79, 137]
n_vals = [1, 2, 2]
l_vals = [0, 0, 1]
kappa_vals = [-1, -1, 1]
p_grid = np.linspace(0, 1000, 100000)

figname = 'comparison_momentum_distribution'
figtitle = 'Comparison between Schrödinger and Dirac momentum distribution functions'
savedir = '../data/' + figname
# %%
# create array of axes in one figure
fig, ax = plt.subplots(3, 3, figsize=(12, 12), constrained_layout=True, gridspec_kw={'width_ratios': [1, 1, 1]})
#fig.suptitle(figtitle, fontsize=22)

for row_index in range(0, 3):
    for col_index in range(0, 3):
        # pick correct set of quantum numbers
        n_val = n_vals[col_index]
        l_val = l_vals[col_index]
        kappa_val = kappa_vals[col_index]
        Z_val = Z_vals[row_index]

        axis = ax[row_index, col_index]

        #manage labeling with axis interface
        if col_index == 2:
            axis.yaxis.set_label_coords(1.0, .5)
            axis.set_ylabel('$Z = {0}$'.format(Z_val), rotation=270, fontsize=22)

        # pick the indexed axis as the current plot
        plt.sca(axis)

        # manage axis labeling
        if row_index == 1 and col_index == 0:
            plt.ylabel(r'$\rho(p)$', fontsize=25)

        if row_index == 2 and col_index == 1:
            plt.xlabel('$p \hspace{5 mm} (a.u.)$', fontsize=25)


        if row_index == 0:
            plt.title('$n = {0}, l = {1}, \kappa = {2}$'.format(n_val, l_val, kappa_val))
        # construct and evaluate Schrödinger momentum distribution function
        dist_nr = calc_momentum_distribution_function(n_val, l_val, Z_val, p_grid)

        # construct and evaluate Dirac momentum distribution function
        dist_rel = calc_momentum_distribution_function_rel(n_val, kappa_val, Z_val, p_grid)

        print(type(dist_rel))

        # rescale axes
        p_grid_rescaled, dist_nr_rescaled, rescale_index = rescale_axes(p_grid, dist_nr, isovalue=0.999)

        dist_rel_rescaled = dist_rel[:rescale_index]

        #plt.plot(p_grid, dist_nr)
        #plt.plot(p_grid, dist_rel)

        plt.plot(p_grid_rescaled, dist_nr_rescaled)
        plt.plot(p_grid_rescaled, dist_rel_rescaled, linestyle='--')

        plt.yscale('log')
        plt.ylim([1e-08, 10])
    
#plt.savefig(fname=savedir, dpi = 300)


plt.show()
# %%
plt.show()
# %%
