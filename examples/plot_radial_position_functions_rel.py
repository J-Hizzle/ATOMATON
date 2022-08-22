# %%
from atomaton.position_space import calc_radial_position_function
from atomaton.position_space_rel import calc_radial_position_function_rel
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
Z_vals = [1, 80, 137]
r_grid = np.linspace(0, 20, 10000)
n_vals = [1, 2, 3]
l_vals = [0, 1, 2]
kappa_vals = [-1, -2, -3]
expression = 'hyp'
# %%
# plot one function for testing purposes 
# define parameters for this loop
#n = 1
#l = 0
#kappa = -1
#Z = 1

#print('n, l, kappa, Z =', n, l, kappa, Z)

#r_grid_rescaled = r_grid

# plot correct functions 
#P_vals_NR = calc_radial_position_function(n, l, Z, r_grid_rescaled)
#P_vals_rel, Q_vals_rel = calc_radial_position_function_rel(n, kappa, Z, r_grid_rescaled, expr=expression)

# write label
#plotlabel_P_NR = r'$P^{NR}_{' + '{0}, {1}'.format(n, l) + '}(r)$'
#plotlabel_P_rel = r'$P^{rel}_{' + '{0}, {1}'.format(n, kappa) + '}(r)$'
#plotlabel_Q_rel = r'$Q^{rel}_{' + '{0}, {1}'.format(n, kappa) + '}(r)$'

#plt.plot(r_grid_rescaled * Z, P_vals_NR/np.sqrt(Z), label=plotlabel_P_NR)
#plt.plot(r_grid_rescaled * Z, P_vals_rel/np.sqrt(Z), label=plotlabel_P_rel)
#plt.plot(r_grid_rescaled * Z, Q_vals_rel/np.sqrt(Z), label=plotlabel_Q_rel)

#plt.ylabel(r'Radial function $(a.u.)$')

#plt.xlabel('$Zr \hspace{5mm} (a.u.)$')

#plt.legend()







# %%
figname = 'Comparison_position_Schr_vs_Dirac_{0}'.format(expression)
figtitle = r'Comparison between Schrödinger and Dirac radial functions for $Z \in \{1, 80, 137\}$ with ' + expression
savedir = '../data/' + figname

# create array of axes in one figure
fig, ax = plt.subplots(3, 3, figsize=(12, 12), constrained_layout=True, gridspec_kw={'width_ratios': [1, 2, 3]})
fig.suptitle(figtitle, fontsize=22)

for row_index in range(0, 3):
    for col_index in range(0, 3):
        # define parameters for this loop
        n = n_vals[col_index]
        l = l_vals[col_index]
        kappa = kappa_vals[col_index]
        Z = Z_vals[row_index]

        #print('n, l, kappa, Z =', n, l, kappa, Z)

        plt.sca(ax[row_index, col_index])

        # manage labeling only the left-most y axes and the bottom-most x axes
        if col_index == 0:
            plt.ylabel(r'Radial function $(a.u.)$')

        if row_index == 2:
            plt.xlabel('$Zr \hspace{5mm} (a.u.)$')

        # define rescaling factor that accounts for distortion at higher quantum numbers
        #rescale = (col_index + 1) * 1/(100 * row_index + 1)

        # rescale r_grid accordingly
        #r_grid_rescaled = r_grid * rescale
        
    	# plot correct functions 
        P_vals_NR = calc_radial_position_function(n, l, Z, r_grid)
        P_vals_rel, Q_vals_rel = calc_radial_position_function_rel(n, kappa, Z, r_grid, expr=expression)

        # write label
        plotlabel_P_NR = r'$P^{NR}_{' + '{0}, {1}'.format(n, l) + '}(r)$'
        plotlabel_P_rel = r'$P^{rel}_{' + '{0}, {1}'.format(n, kappa) + '}(r)$'
        plotlabel_Q_rel = r'$Q^{rel}_{' + '{0}, {1}'.format(n, kappa) + '}(r)$'
        
        plt.plot(r_grid * Z, P_vals_NR/np.sqrt(Z), label=plotlabel_P_NR)
        plt.plot(r_grid * Z, P_vals_rel/np.sqrt(Z), label=plotlabel_P_rel)
        plt.plot(r_grid * Z, Q_vals_rel/np.sqrt(Z), label=plotlabel_Q_rel)

        plt.legend('yop')

# postprocess figure as a whole
# manage legend for the entire figure
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
handles, labels = lines_labels[0][0], lines_labels[0][1]
fig.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5))

#plt.savefig(fname=savedir, dpi = 300)

# %%
