# %%
import sympy as smp
# %%
'''
=================
= Preliminaries =
=================
'''
# define primary variables, parameters and constants and specify their properties
r = smp.symbols('r', real=True, positive=True)                      # radial position
n, Z = smp.symbols('n, Z', integer=True, positive=True)             # principal quantum number, nuclear charge number
kappa = smp.symbols('kappa', integer=True, nonzero=True)            # relativistic symmetry quantum number
alpha = smp.symbols('alpha', positive=True)                         # inverse speed of light

# define secondary variables, parameters and constants 
gamma = smp.sqrt(kappa**2 - Z**2 * alpha**2)                        # apparent angular quantum number
n_r = n - smp.Abs(kappa)                                            # radial quantum number
N = smp.sqrt(n_r**2 + 2 * n_r * gamma + kappa**2)                   # apparent principal quantum number
C_P = smp.sqrt(1 + (n_r + gamma)/N)                                 # prefactor of large component
C_Q = smp.sqrt(1 - (n_r + gamma)/N)                                 # prefactor of small component
norm_lag = smp.sqrt(Z * smp.factorial(n_r)/(2 * N**2 * (N - kappa) * smp.gamma(2 * gamma + n_r + 1)))   # normalization constant
norm_hyp = smp.sqrt(Z * smp.factorial(n_r + 2 * gamma)/(2 * N**2 * (N - kappa) * smp.factorial(2 * gamma)**2 * smp.factorial(n_r)))
x = 2 * Z/N * r                                                  # rescaled radial variable

# define radial functions P and Q in Andrae's representation with Laguerre polynomials as a function of the radius r in atomic units
P_lag = norm_lag * C_P * x**gamma * smp.exp(-x/2) * ((N - kappa) * smp.assoc_laguerre(n_r, 2 * gamma, x) - n_r * smp.assoc_laguerre(n_r - 1, 2 * gamma, x))
Q_lag = -norm_lag * C_Q * x**gamma * smp.exp(-x/2) * ((N - kappa) * smp.assoc_laguerre(n_r, 2 * gamma, x) + n_r * smp.assoc_laguerre(n_r - 1, 2 * gamma, x))

# define special case for n_r = 0
P_lag_0 = norm_lag * C_P * x**gamma * smp.exp(-x/2) * (N - kappa) * smp.assoc_laguerre(n_r, 2 * gamma, x)
Q_lag_0 = -norm_lag * C_Q * x**gamma * smp.exp(-x/2) * (N - kappa) * smp.assoc_laguerre(n_r, 2 * gamma, x)


# define radial functions P and Q in Andrae's representation with confluent hypergeometric functions as a function of the radius r in atomic units
P_hyp = norm_hyp * C_P * x**gamma * smp.exp(-x/2) * ((N - kappa) * smp.hyper([-n_r], [2 * gamma + 1], x) - n_r * smp.hyper([-n_r + 1], [2 * gamma + 1], x))
Q_hyp = -norm_hyp * C_Q * x**gamma * smp.exp(-x/2) * ((N - kappa) * smp.hyper([-n_r], [2 * gamma + 1], x) + n_r * smp.hyper([-n_r + 1], [2 * gamma + 1], x))

# define special case for n_r = 0
#P_hyp_0 = norm_lag * C_P * x**gamma * smp.exp(-x/2) * (N - kappa) * smp.hyper([-n_r], [2 * gamma + 1], x)
#Q_hyp_0 = -norm_lag * C_Q * x**gamma * smp.exp(-x/2) * (N - kappa) * smp.hyper([-n_r], [2 * gamma + 1], x)
# %%
'''
=============
= Functions =
=============
'''
def calc_radial_position_function_rel(n_val, kappa_val, Z_val, r_grid, alpha_val=1/137.035999037, expr='lag'):
    '''
    Construct radial position position function P_{n, l}(r) in a.u. and evaluate at radii given by r_grid.
    '''
    if expr == 'lag':
        # intercept error for negative degree of laguerre polynomials
        if n_val - (kappa_val**2)**(1/2) == 0:
    
            # subsitute in given parameters
            P_nl = P_lag_0.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])
            Q_nl = Q_lag_0.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])

        else:
            # subsitute in given parameters
            P_nl = P_lag.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])
            Q_nl = Q_lag.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])

    elif expr == 'hyp':
        P_nl = P_hyp.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])
        Q_nl = Q_hyp.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])
    
    P_nl = smp.simplify(P_nl)
    Q_nl = smp.simplify(Q_nl)
    
    # lambdify sympy expression
    P_func = smp.lambdify(r, P_nl)
    Q_func = smp.lambdify(r, Q_nl)

    return P_func(r_grid), Q_func(r_grid)
#%%
#def calc_r_k_expval(n_val, l_val, Z_val, k_val):
    '''
    Compute expectation value of r^k as integral
    '''
    # substitute in given parameters
    #P_nl = P.subs([(n, n_val), (l, l_val), (Z, Z_val)])

    # calculate integral
    #r_k_expval = smp.integrate(P_nl**2 * r**k_val, (r, 0, smp.oo))

    #return r_k_expval