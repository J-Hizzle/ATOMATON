# %%
import sympy as smp
# %%
'''
=================
= Preliminaries =
=================
'''
# define all relevant variables and specify their properties
p = smp.symbols('p', real=True, positive=True)                      # radial momentum
n, Z = smp.symbols('n, Z', integer=True, positive=True)             # principal quantum number, nuclear charge number
l = smp.symbols('l', integer=True, nonnegative=True)                # angular momentum quantum number
k = smp.symbols('k', integer=True)                                  # exponent of p expectation value
t = smp.symbols('t', real=True)                                     # transformed radial momentum variable for integration

# define radial momentum function G as a function of the radial momentum p in atomic units with phase factor omitted
G = n * 2**(2*l + 2) * smp.factorial(l) * smp.sqrt(2 * smp.factorial(n - l - 1)/(Z * smp.pi * smp.factorial(n + l))) \
    * (p * n/Z)**(l + 1)/((p * n/Z)**2 + 1)**(l + 2) * smp.gegenbauer(n - l - 1, l + 1, ((p * n/Z)**2 - 1)/((p * n/Z)**2 + 1))

# define radial part F as quotient of radial function G and radial momentum p
F = G/p

# define integrand of expectation values 
i = smp.gegenbauer(n - l - 1, l + 1, t)**2 * (1 - t)**(l + 3/2 - 1/2 * k) * (1 + t)**(l + 1/2 + 1/2 * k)
# %%
'''
=============
= Functions =
=============
'''
def calc_real_radial_momentum_function(n_val, l_val, Z_val, p_grid):
    '''
    Construct radial momentum function G_{n, l}(p) with in a.u. and evaluate at radial momenta given by p_grid.
    '''
    # subsitute parameters into sympy expression
    G_nl = G.subs([(n, n_val), (l, l_val), (Z, Z_val)])

    # lambdify sympy expression
    radial_mom = smp.lambdify(p, G_nl)

    return radial_mom(p_grid)
# %%
def calc_p_k_expval(n_val, l_val, Z_val, k_val):
    '''
    Compute expectation value of p^k as integral
    '''
    # substitute in given parameters
    G_nl = G.subs([(n, n_val), (l, l_val), (Z, Z_val)])

    # calculate integral
    p_k_expval = smp.integrate(G_nl**2 * p**k_val, (p, 0, smp.oo))

    return p_k_expval
# %%
def calc_expval_integrand(n_val, l_val, k_val, t_grid):
    '''
    Compute integrand of the expectation value calculation to showcase its properties and evaluate at given t_grid
    '''
    # substitute in given parameters
    i_knl = i.subs([(n, n_val), (l, l_val), (k, k_val)])

    # lambdify sympy expression
    integrand = smp.lambdify(t, i_knl)

    return integrand(t_grid)
# %%