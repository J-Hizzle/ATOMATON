# %%
import sympy as smp
# %%
'''
=================
= Preliminaries =
=================
'''
# define all relevant variables and specify their properties
r = smp.symbols('r', real=True, positive=True)                      # radial position
n, Z = smp.symbols('n, Z', integer=True, positive=True)             # principal quantum number, nuclear charge number
l = smp.symbols('l', integer=True, nonnegative=True)                # angular momentum quantum number

# define radial part R as a function of the radius r in atomic units
R = smp.sqrt((2 * Z/n)**3 * smp.factorial(n-l-1)/(2*n*(smp.factorial(n+l)))) \
    * smp.exp(-Z * r/n) * (2 * Z * r/n)**l * smp.assoc_laguerre(n-l-1,2*l+1,(2 * Z * r/n))

# define radial position function P as product of radial part R and radius r
P = R * r
# %%
'''
=============
= Functions =
=============
'''
def calc_radial_position_function(n_val, l_val, Z_val, r_grid):
    '''
    Construct radial position position function P_{n, l}(r) in a.u. and evaluate at radii given by r_grid.
    '''
    # subsitute substitute in given parameters
    P_nl = P.subs([(n, n_val), (l, l_val), (Z, Z_val)])

    # lambdify sympy expression
    radial_dist = smp.lambdify(r, P_nl)

    return radial_dist(r_grid)
#%%
def calc_radial_distribution_function(n_val, l_val, Z_val, r_grid):
    '''
    Construct and evaluate momentum distribution function.
    '''
    # subsitute parameters into sympy expression
    P_nl = P.subs([(n, n_val), (l, l_val), (Z, Z_val)])

    # lambdify sympy expression
    P_func = smp.lambdify(r, P_nl)

    return P_func(r_grid)**2
# %%
def calc_r_k_expval(n_val, l_val, Z_val, k_val):
    '''
    Compute expectation value of r^k as integral
    '''
    # substitute in given parameters
    P_nl = P.subs([(n, n_val), (l, l_val), (Z, Z_val)])

    # calculate integral
    r_k_expval = smp.integrate(P_nl**2 * r**k_val, (r, 0, smp.oo))

    return r_k_expval