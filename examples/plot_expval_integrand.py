'''
=====================
plot_expval_integrand
=====================

Showcase how integrand i_knl of expectation value integral looks like for different k

1) Evaluate functions for certain combinations of n, k

2) Plot solutions with matplotlib
'''
# %%
from turtle import position
from matplotlib.lines import _LineStyle
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
n_val = 1
l_val = 0
#k_vals = np.arange(-3, 6, 1, dtype=int)
k_vals = np.flip(np.array([-3, -2, -1, 3, 4, 5]))
t_grid = np.linspace(-1+1e-10, 1, 1000000)

# specify names and save path
figname = 'expval_integrands'
figtitle = r'Integrands $i_{k,n,l}(t)$ vs integration variable $t$ for $n = 1, l = 0$'
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
exp_vals = np.zeros((len(k_vals), len(t_grid)))

# calculate expectation values
for jj, k_val in enumerate(k_vals):
    exp_vals[jj, :] = calc_expval_integrand(n_val, l_val, k_val, t_grid)
    print('calc state =', jj)
# %%
'''
======
= 3) =
======
'''

plt.figure(figsize=(12, 12))
plt.ylabel(r'$i_{k,n,l}(t) \hspace{5mm} (a.u.)$')
plt.xlabel('$t \hspace{5mm} (a.u.)$')

# plot integrand values i_knl against integration variable t for multiple k_vals
for jj, k_val in enumerate(k_vals):
    plt.plot(t_grid, exp_vals[jj, :], label='$k = {}$'.format(k_val))

# calculate auxiliary functions
f_plu = f_plus(t_grid)
f_min = f_minus(t_grid)

# plot function 1/(x - 1) and 1/(x + 1) as reference
plt.plot(t_grid, f_plu, '--', label='$1/(1 + t)$', linewidth=1.5, color='black')
plt.plot(t_grid, f_min, '--', label='$1/t - t)$', linewidth=1.5, color='black')

plt.xticks([-1.0, -0.5, 0.0, 0.5, 1.0], [-1.0, -0.5, 0.0, 0.5, 1.0])
plt.xlim(-1.2, 1.2)

# define list of ticks for the y axis
#y_tick_list = np.asarray([10**x for x in range(-5, 13)])

#plt.yticks(y_tick_list, y_tick_list)
plt.ylim(-20, 100)
#plt.yscale('log')

plt.title(figtitle)
plt.legend(loc=10)

#plt.tick_params(which='minor', width=0)

# postprocess figure as a whole
# manage legend for the entire figure
#lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
#handles, labels = lines_labels[0][0], lines_labels[0][1]
#fig.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5))

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
