# %%
import sympy as smp
import numpy as np
# %%
'''
=================
= Preliminaries =
=================
'''
# define primary variables, parameters and constants and specify their properties
p = smp.symbols('p', real=True, positive=True)                      # radial momentum 
n, Z = smp.symbols('n, Z', integer=True, positive=True)             # principal quantum number, nuclear charge number
kappa = smp.symbols('kappa', integer=True, nonzero=True)            # relativistic symmetry quantum number
alpha = smp.symbols('alpha', positive=True)                         # inverse speed of light
k, j = smp.symbols('k, j', integer=True)                            # index number for summation

# define secondary variables, parameters and constants 
a_s = smp.sign(kappa)                                               # negative sign factor
l = a_s * (kappa - 1/2 * (1 + a_s))                                 # angular quantum number of upper component
l_prime = l + a_s                                                   # angular quantum number of lower component
gamma = smp.sqrt(kappa**2 - Z**2 * alpha**2)                        # apparent angular quantum number
n_r = n - smp.Abs(kappa)                                            # radial quantum number
N = smp.sqrt(n_r**2 + 2 * n_r * gamma + kappa**2)                   # apparent principal quantum number
C_plus = smp.sqrt(1 + (n_r + gamma)/N)                              # prefactor of large component
C_minus = smp.sqrt(1 - (n_r + gamma)/N)                             # prefactor of small component
M_norm_plus = smp.I**l * 2**(gamma - l + 2)/smp.gamma(l + 3/2) * smp.sqrt((N**2 * smp.gamma(n_r + 2 * gamma + 1)))/(Z**3 * (N - kappa) * smp.gamma(2 * gamma + 1) * smp.factorial(n_r))                     # normalization constant for the upper component
M_norm_minus = smp.I**l_prime * 2**(gamma - l_prime + 2)/smp.gamma(l_prime + 3/2) * smp.sqrt((N**2 * smp.gamma(n_r + 2 * gamma + 1))/(Z**3 * (N - kappa) * smp.gamma(2 * gamma + 1) * smp.factorial(n_r)))  # normalization constant for the lower component
y = N/Z * p                                                         # rescaled radial variable

# define subfunctions contained in radial momentum functions
S_1_plus = smp.rf(-n_r, kappa)/(smp.rf(2 * gamma + 1, k) * 2**k * smp.factorial(k)) * y**l/(y**2 + 1)**(1/2 * (l + gamma + k + 2)) * smp.hyper([(l + gamma + k + 2)/2, (l - gamma - k)/2], [(2 * l + 3)/2], y**2/(y**2 + 1))  
S_2_plus = smp.rf(-n_r + 1, kappa)/(smp.rf(2 * gamma + 1, k) * 2**k * smp.factorial(k)) * y**l/(y**2 + 1)**(1/2 * (l + gamma + k + 2)) * smp.hyper([(l + gamma + k + 2)/2, (l - gamma - k)/2], [(2 * l + 3)/2], y**2/(y**2 + 1))
S_1_minus = smp.rf(-n_r, kappa)/(smp.rf(2 * gamma + 1, k) * 2**k * smp.factorial(k)) * y**l_prime/(y**2 + 1)**(1/2 * (l_prime + gamma + k + 2)) * smp.hyper([(l_prime + gamma + k + 2)/2, (l_prime - gamma - k)/2], [(2 * l_prime + 3)/2], y**2/(y**2 + 1))  
S_2_minus = smp.rf(-n_r + 1, kappa)/(smp.rf(2 * gamma + 1, k) * 2**k * smp.factorial(k)) * y**l_prime/(y**2 + 1)**(1/2 * (l_prime + gamma + k + 2)) * smp.hyper([(l_prime + gamma + k + 2)/2, (l_prime - gamma - k)/2], [(2 * l_prime + 3)/2], y**2/(y**2 + 1))

G_1_plus = smp.summation(S_1_plus, (k, 0, n_r))
G_2_plus = smp.summation(S_2_plus, (k, 0, n_r - 1))
G_1_minus = smp.summation(S_1_minus, (k, 0, n_r))
G_2_minus = smp.summation(S_2_minus, (k, 0, n_r - 1))

# define radial momentum functions G and H obtained from Fourier transformation of position space spinor
G = M_norm_plus * C_plus * ((N - kappa) * G_1_plus - n_r * G_2_plus)
H = M_norm_minus * C_minus * ((N - kappa) * G_1_minus + n_r * G_2_minus)
# %%
def calc_momentum_distribution_function_rel(n_val, kappa_val, Z_val, p_grid, alpha_val=1/137.035999037):
    '''
    Construct and evaluate radial momentum distribution function.
    '''
    G_nk = G.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])
    H_nk = H.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])

    #G_nk = smp.simplify(G_nk)
    #H_nk = smp.simplify(H_nk)

    dist = G_nk**2 + H_nk**2

    dist_func = []

    # evaluate sympy expressions numerically with looping over p_grid

    for ii, p_val in enumerate(p_grid):
        dist_func.append(dist.evalf(subs={p: p_val}))
        print('ii =', ii)
        print('dist_func =', dist_func[ii])

    dist_func = np.asarray(dist_func)

    # lambdify sympy expressions
    #G_func = smp.lambdify(p, G_nk)
    #H_func = smp.lambdify(p, H_nk)

    #print('type(G_func) =', type(G_func))

    print('type(dist_func) =', type(dist_func))

    return dist_func
# %%
if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt

    summation = G_dunno**2

    sumsub = summation.subs([(n, 1), (kappa, -1), (Z, 1), (alpha, 1/137.035999037)])

    #sumsub = smp.simplify(sumsub)

    p_grid = np.linspace(0, 20, 100)

    sum_vals = []

    for ii, p_val in enumerate(p_grid):
        sum_vals.append(sumsub.evalf(subs={p: p_val}))
        print('ii =', ii)
        print('type(val) =', type(sum_vals[ii]))
        print('val =', sum_vals[ii])

    G_vals = np.asarray(sum_vals)

    plt.plot(p_grid, sum_vals)
    # %%
    plt.show()

    # %%
    dist_func = calc_momentum_distribution_function_rel(1, -1, 1, p_grid)
# %%