#%%
import sympy as smp
import numpy as np
import matplotlib.pyplot as plt
from atomaton.momentum_space import calc_p_k_expval
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
'''
=========
p_expvals
=========

Get lists of expectation values p^k for k = -2, -1, 0, 1, 2, 3, 4 from integration

1) Evaluate integrals for certain combinations of n, l
2a) Evaluate variance of r^1 expectation value

3) Plot solutions with matplotlib
'''
# %%
'''
=============
preliminaries
=============
'''
# define array over principal quantum numbers and set all other parameters constant
n_vals = np.arange(1, 10, 1, dtype=int)
l_vals = np.arange(0, 4, dtype=int)
Z_val = 1
k_vals = np.arange(-2, 5, 1, dtype=int)

# specify names and save path
figname = 'p_exp_vals_Z{}_nocap'.format(Z_val)
figtitle = r'Radial momentum expectation values $\langle p^k \rangle$ vs. principal quantum numbers for $Z = {}$'.format(Z_val)
savedir = '../data/' + figname
# %%
'''
======
= 1) =
======
'''
# initialize array for results
exp_vals = np.zeros((len(l_vals), len(k_vals), len(n_vals)))

# calculate expectation values
for kk, l_val in enumerate(l_vals):
    for jj, k_val in enumerate(k_vals):
        for ii, n_val in enumerate(n_vals):
            if l_val < n_val:
                exp_vals[kk, jj, ii] = calc_p_k_expval(n_val, l_val, Z_val, k_val)
            else: 
                exp_vals[kk, jj, ii] = np.NaN
            print('calc state =', kk, jj, ii)
# %%
'''
======
= 3) =
======
'''
# initiate number of desired subplots in each direction
row_number = 2
col_number = len(l_vals) - row_number
l_index = 0

# create array of axes in one figure
fig, ax = plt.subplots(row_number, col_number, figsize=(12, 12), constrained_layout=True, gridspec_kw={'width_ratios': [1, 1]})
#fig.suptitle(figtitle, fontsize=22)

# loop over subplots in figure and fill them with data 
for row_index in range(0, row_number):
    for col_index in range(0, col_number):
        # determine which value of l is the current one
        l_val = l_vals[l_index]

        # pick one of the subplots as the current plot
        plt.sca(ax[row_index, col_index])

        # manage labeling only the left-most y axes and the bottom-most x axes
        #if col_index == 0:
        #    plt.ylabel(r'$\langle p^k \rangle \hspace{5mm} (a.u.)$')
            #plt.xlim(0, 4)
        #elif col_index == col_number - 1:
            #plt.xlim(0, 2)

        if row_index == row_number - 1:
            plt.xlabel('$n$')

        # plot expectation values against principal quantum numbers
        for jj, k_val in enumerate(k_vals):
            plt.scatter(n_vals, exp_vals[l_val, jj, :], label='$k = {}$'.format(k_val))
            
        plt.xticks(n_vals, n_vals)
        plt.xlim(n_vals[0] - 1, n_vals[-1] + 1)
        
        # define list of ticks for the y axis
        y_tick_list = np.asarray([10**x for x in range(-5, 13)])

        plt.yticks(y_tick_list, y_tick_list)
        plt.ylim(1e-5, 1e12)
        plt.yscale('log')

        plt.title('$l = {}$'.format(l_val))
        plt.tick_params(which='minor', width=0)
        # index over angular momentum quantum numbers
        l_index += 1

# postprocess figure as a whole
# manage legend for the entire figure
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
handles, labels = lines_labels[0][0], lines_labels[0][1]
fig.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5))

# hide every second label of the x ticks
for ax in fig.axes:
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)

# hide every second label of the y ticks
#for ax in fig.axes:
#    for label in ax.yaxis.get_ticklabels()[::2]:
#        label.set_visible(False)

# label x and y axes for entire figure
fig.supylabel(r'$\langle p^k \rangle \hspace{5mm} (a.u.)$', fontsize=25)
fig.supxlabel('$n$', fontsize=25)

    
plt.savefig(fname=savedir, dpi = 300)
plt.show()
# %%
