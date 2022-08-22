# %%
import numpy as np
import sympy as smp
from sympy import gegenbauer
# %%
# define all relevant variables and specify their properties
p = smp.symbols('p', real=True, positive=True)                      # radial momentum
n, l, Z = smp.symbols('n, l, Z', integer=True, positive=True)       # principal quantum number, angular momentum quantum number, nuclear charge number

# define radial momentum distribution function G as a function of the radial momentum p in atomic units
G = -(-smp.I)**l * n * 2**(2*l + 2) * smp.factorial(l) * smp.sqrt(2 * smp.factorial(n - l - 1)/(Z * smp.pi * smp.factorial(n + l))) \
    * (p * n/Z)**(l + 1)/((p * n/Z)**2 + 1)**(l + 2) * smp.gegenbauer(n - l - 1, l + 1, ((p * n/Z)**2 - 1)/((p * n/Z)**2 + 1))

# define radial part F as quotient of radial function G and radial momentum p
F = G/p

# %%
def calc_real_radial_momentum_function(n_val, l_val, Z_val, p_grid):
    '''
    Construct radial momentum distribution function G_{n, l}(p) with in a.u. and evaluate at radial momenta given by p_grid.
    '''
    # subsitute parameters into sympy expression
    G_nl = smp.re(G.subs([(n, n_val), (l, l_val), (Z, Z_val)]))
    print('G_nl =', G_nl)

    G_nl_alt = G.subs([(n, n_val), (l, l_val), (Z, Z_val)])
    print('G_nl_alt =', G_nl_alt)


    # lambdify sympy expression
    radial_dist = smp.lambdify(p, G_nl)

    return radial_dist(p_grid)


# %%
def radial_momentum_function_alt(n, l, Z, y_grid):
    '''
    [deprecated]
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
