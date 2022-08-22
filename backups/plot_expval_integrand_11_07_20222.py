'''
=====================
plot_expval_integrand
=====================

Showcase how integrand i_knl of expectation value integral looks like for different k

1) Evaluate functions for certain combinations of n, k

2) Plot solutions with matplotlib
'''
# %%
import numpy as np
import matplotlib.pyplot as plt
from atomaton.momentum_space import calc_expval_integrand
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title

# %%
'''
=============
preliminaries
=============
'''
# define array over principal quantum numbers and set all other parameters constant
n_vals = np.arange(1, 5, 1, dtype=int)
l_val = 0
#k_vals = np.arange(-3, 6, 1, dtype=int)
k_vals = np.array([-3, -2, -1, 3, 4, 5])

t_grid = np.linspace(-1+1e-10, 1, 1000000)

# specify names and save path
figname = 'expval_integrands'
figtitle = r'Integrands $i_{k,n,l}(t)$ vs integration variable $t$ appearing in $\langle p^k \rangle$'
savedir = '../data/' + figname


# %%
# define auxiliary functions 
def f_minus(t_grid):
    f = np.zeros(len(t_grid))

    for ii, t_val in enumerate(t_grid):
        if 1 - np.abs(t_val) > 1e-10:
            f[ii] = 1/(1 - t_val)
        else:
            f[ii] = np.NAN
    
    return f

def f_plus(t_grid):
    f = np.zeros(len(t_grid))

    for ii, t_val in enumerate(t_grid):
        if 1 - np.abs(t_val) > 1e-10:
            f[ii] = 1/(1 + t_val)
        else:
            f[ii] = np.NAN
    
    return f
# %%
'''
======
= 1) =
======
'''
# initialize array for results
exp_vals = np.zeros((len(n_vals), len(k_vals), len(t_grid)))

# calculate expectation values
for kk, n_val in enumerate(n_vals):
    for jj, k_val in enumerate(k_vals):
        exp_vals[kk, jj, :] = calc_expval_integrand(n_val, l_val, k_val, t_grid)
        print('calc state =', kk, jj)
# %%
'''
======
= 3) =
======
'''
# initiate number of desired subplots in each direction
row_number = 2
col_number = len(n_vals) - row_number
n_index = 0

# create array of axes in one figure
fig, ax = plt.subplots(row_number, col_number, figsize=(12, 12), constrained_layout=True, gridspec_kw={'width_ratios': [1, 1]})
fig.suptitle(figtitle, fontsize=22)

# loop over subplots in figure and fill them with data 
for row_index in range(0, row_number):
    for col_index in range(0, col_number):
        # determine which value of l is the current one
        n_val = n_vals[n_index]

        # pick one of the subplots as the current plot
        plt.sca(ax[row_index, col_index])

        # manage labeling only the left-most y axes and the bottom-most x axes
        if col_index == 0:
            plt.ylabel(r'$i_{k,n,l}(t) \hspace{5mm} (a.u.)$')
            #plt.xlim(0, 4)
        #elif col_index == col_number - 1:
            #plt.xlim(0, 2)

        if row_index == row_number - 1:
            plt.xlabel('$t \hspace{5mm} (a.u.)$')

        # plot integrand values i_knl against integration variable t for multiple k_vals
        for jj, k_val in enumerate(k_vals):
            plt.plot(t_grid, exp_vals[n_index, jj, :], label='$k = {}$'.format(k_val))

        # calculate auxiliary functions
        f_plu = f_plus(t_grid)
        f_min = f_minus(t_grid)

        # plot function 1/(x - 1) and 1/(x + 1) as reference
        plt.plot(t_grid, f_plu, 'r--', label='$1/(1 + t)$')
        plt.plot(t_grid, f_min, 'r--', label='$1/t - t)$')

        plt.xticks([-1.0, -0.5, 0.0, 0.5, 1.0], [-1.0, -0.5, 0.0, 0.5, 1.0])
        plt.xlim(-1.2, 1.2)
        
        # define list of ticks for the y axis
        #y_tick_list = np.asarray([10**x for x in range(-5, 13)])

        #plt.yticks(y_tick_list, y_tick_list)
        plt.ylim(1e-10, 5e2)
        #plt.yscale('log')

        plt.title('$n = {}$'.format(n_val))
        plt.tick_params(which='minor', width=0)
        # index over angular momentum quantum numbers
        n_index += 1

# postprocess figure as a whole
# manage legend for the entire figure
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
handles, labels = lines_labels[0][0], lines_labels[0][1]
fig.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5))

# hide every second label of the x ticks
#for ax in fig.axes:
#    for label in ax.xaxis.get_ticklabels()[::2]:
#        label.set_visible(False)

# hide every second label of the y ticks
#for ax in fig.axes:
#    for label in ax.yaxis.get_ticklabels()[::2]:
#        label.set_visible(False)
    
plt.savefig(fname=savedir, dpi = 300)
plt.show()
# %%
