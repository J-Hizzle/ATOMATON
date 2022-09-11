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
k, j, i, h = smp.symbols('k, j, i, h', integer=True)                            # index number for summation

#s econdary variables
a_s = -smp.sign(kappa)    
l = -(a_s * kappa + 1/2 + 1/2 * a_s)
l_prime = l + a_s  
n_r = n - smp.Abs(kappa)
a = alpha * Z
gamma = smp.sqrt(kappa**2 - a**2) 
N = smp.sqrt(n_r**2 + 2 * n_r * gamma + kappa**2)                   # apparent principal quantum number
lamb = Z/N
W = smp.sqrt(1 - lamb**2)
#1/alpha**2 * (1/smp.sqrt(1 + (a/(n_r + gamma)**2)) - 1)
C_n_r = -2**(gamma - 1/2) * lamb**(gamma + 3/2) * smp.sqrt(smp.gamma(2 * gamma + n_r + 1)/(smp.factorial(n_r) * a * (a - lamb * kappa)))
# %%
# define radial momentum summands
S_A_l = smp.rf(-n_r + 1, k) * smp.gamma(l + k + gamma + 2) * (2 * lamb)**k/(smp.gamma(k + 2 * gamma + 1) * smp.gamma(k + 1) * (lamb**2 + p**2)**(1/2 * (k + gamma + l + 2))) * smp.hyper([(k + gamma + l + 2)/2, (l - k - gamma)/2], [(2 * l + 3)/2], p**2/(lamb**2 + p**2))
S_B_l = smp.rf(-n_r, k) * smp.gamma(l + k + gamma + 2) * (2 * lamb)**k/(smp.gamma(k + 2 * gamma + 1) * smp.gamma(k + 1) * (lamb**2 + p**2)**(1/2 * (k + gamma + l + 2))) * smp.hyper([(k + gamma + l + 2)/2, (l - k - gamma)/2], [(2 * l + 3)/2], p**2/(lamb**2 + p**2))
S_A_l_prime = smp.rf(-n_r + 1, k) * smp.gamma(l_prime + k + gamma + 2) * (2 * lamb)**k/(smp.gamma(k + 2 * gamma + 1) * smp.gamma(k + 1) * (lamb**2 + p**2)**(1/2 * (k + gamma + l_prime + 2))) * smp.hyper([(k + gamma + l_prime + 2)/2, (l_prime - k - gamma)/2], [(2 * l_prime + 3)/2], p**2/(lamb**2 + p**2))
S_B_l_prime = smp.rf(-n_r, k) * smp.gamma(l_prime + k + gamma + 2) * (2 * lamb)**k/(smp.gamma(k + 2 * gamma + 1) * smp.gamma(k + 1) * (lamb**2 + p**2)**(1/2 * (k + gamma + l_prime + 2))) * smp.hyper([(k + gamma + l_prime + 2)/2, (l_prime - k - gamma)/2], [(2 * l_prime + 3)/2], p**2/(lamb**2 + p**2))
# %%
# define radial momentum function parts
A_l = smp.gamma(1/2) * p**l/(2**(l + 1) * smp.gamma(l + 3/2)) * smp.summation(S_A_l, (k, 0, n_r - 1))
B_l = smp.gamma(1/2) * p**l/(2**(l + 1) * smp.gamma(l + 3/2)) * smp.summation(S_B_l, (k, 0, n_r))
A_l_prime = smp.gamma(1/2) * p**l_prime/(2**(l_prime + 1) * smp.gamma(l_prime + 3/2)) * smp.summation(S_A_l_prime, (k, 0, n_r - 1))
B_l_prime = smp.gamma(1/2) * p**l_prime/(2**(l_prime + 1) * smp.gamma(l_prime + 3/2)) * smp.summation(S_B_l_prime, (k, 0, n_r))
# %%
# define radial momentum functions
N = (2/smp.pi)**(1/2) * (1 + W)**(1/2) * C_n_r * (-n_r * A_l - (kappa - a/lamb) * B_l)
#-smp.I**(-l) * (2/smp.pi)**(1/2) * (1 + W)**(1/2) * C_n_r * (-n_r * A_l - (kappa - a/lamb) * B_l)
M = (2/smp.pi)**(1/2) * (1 - W)**(1/2) * C_n_r * (n_r * A_l_prime - (kappa - a/lamb) * B_l_prime)
#smp.I**(-l_prime) * (2/smp.pi)**(1/2) * (1 - W)**(1/2) * C_n_r * (n_r * A_l_prime - (kappa - a/lamb) * B_l_prime)



# %%
def calc_momentum_distribution_function_rel(n_val, kappa_val, Z_val, p_grid, alpha_val=1/137.035999037):
    '''
    Construct and evaluate radial momentum distribution function.
    '''
    #if n_val - np.abs(kappa_val) == 0:
    #    G_nk = G_alt.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])
    #    H_nk = H_alt.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])

    #else:
    G_nk = N.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])
    H_nk = M.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])

    print(type(G_nk))

    #G_nk = smp.simplify(G_nk)
    #H_nk = smp.simplify(H_nk)

    dist = (G_nk * smp.conjugate(G_nk) + H_nk * smp.conjugate(H_nk)) * p**2

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

    n_val = 1
    kappa_val = -1
    Z_val = 1
    alpha_val = 1/137.035999037
    
    G_nk = N.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])
    H_nk = M.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)])

    #sumsub = smp.simplify(sumsub)

    sumsub = (G_nk * smp.conjugate(G_nk) + H_nk * smp.conjugate(H_nk)) * p**2

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
    smp.evalf(lamb.subs([(n, n_val), (kappa, kappa_val), (Z, Z_val), (alpha, alpha_val)]))