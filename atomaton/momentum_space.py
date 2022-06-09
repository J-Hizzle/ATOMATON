# %%
from scipy.special import gegenbauer, factorial
import numpy as np
# %%
def radial_function(n, l, Z, y_grid):
    '''
    Construct radial function Z_{n, l}(y) with y = np/Z in a.u. and evaluate at position y.
    '''
    # construct gegenbauer polynomial C for given n and l
    C_vals = gegenbauer(n - l - 1, l + 1)((y_grid**2 - 1)/(y_grid**2 + 1))

    # calculate prefactor N
    N = n * 2**(2 * l + 2) * factorial(l) * np.sqrt((2 * factorial(n - l - 1))/(Z * np.pi * factorial(n + l)))
    #print('N =', N)

    # evaluate polynomial part P at positions y
    P_vals = y_grid**(l + 1)/(y_grid**2 + 1)**(l + 2)
    #print('P_vals =', P_vals)

    # put all three parts together 
    Z_vals = N * C_vals * P_vals

    return Z_vals
# %%
def momental(n, l, m, Z, )