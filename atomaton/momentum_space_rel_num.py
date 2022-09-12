# %%
import sympy as smp
import numpy as np
import scipy as scp
import scipy.special as spc
# %%
def calc_momentum_distribution_function_rel(n_val, kappa_val, Z_val, p_grid, alpha_val=1/137.035999037):
    n = n_val
    kappa = kappa_val
    Z = Z_val
    alpha = alpha_val
    
    # calculate secondary variables and constants
    a_s = -np.sign(kappa)
    l = -(a_s * kappa + 1/2 + a_s/2)
    l_prime = l + a_s
    gamma = np.sqrt(kappa**2 - Z**2 * alpha**2)
    n_r = n - np.abs(kappa)
    N = np.sqrt(n_r**2 + 2 * n_r * gamma + kappa**2)
    C_plus = np.sqrt(1 + (n_r + gamma)/N)
    C_minus = np.sqrt(1 - (n_r + gamma)/N)
    M_norm_noimag = 2**(gamma - l - 1)/spc.gamma(l + 3/2) * np.sqrt(spc.gamma(n_r + 2 * gamma + 1)/(Z * (N - kappa) * (spc.gamma(2 * gamma + 1))**2 * spc.factorial(n_r)))
    R_norm_noimag = 2**(gamma - l_prime - 1)/spc.gamma(l_prime + 3/2) * np.sqrt(spc.gamma(n_r + 2 * gamma + 1)/(Z * (N - kappa) * (spc.gamma(2 * gamma + 1))**2 * spc.factorial(n_r)))
    #y_grid = Z/N * p_grid
    y_grid = N/Z * p_grid

    G_1 = np.zeros_like(y_grid)
    G_2 = np.zeros_like(y_grid)
    F_1 = np.zeros_like(y_grid)
    F_2 = np.zeros_like(y_grid)

    for k in range(0, n_r + 1):
        G_1 += spc.gamma(l + gamma + k + 2) * spc.poch(-n_r, k)/spc.poch(2 * gamma + 1, k) * 2**k/spc.factorial(k) * y_grid**(l + 1)/(y_grid**2 + 1)**((l + gamma + k +2)/2) \
            * spc.hyp2f1((l + gamma + k + 2)/2, (l - gamma - k)/2, (2 * l + 3)/2, y_grid**2/(y_grid**2 + 1))
        F_1 += spc.gamma(l_prime + gamma + k + 2) * spc.poch(-n_r, k)/spc.poch(2 * gamma + 1, k) * 2**k/spc.factorial(k) * y_grid**(l_prime + 1)/(y_grid**2 + 1)**((l_prime + gamma + k +2)/2) \
            * spc.hyp2f1((l_prime + gamma + k + 2)/2, (l_prime - gamma - k)/2, (2 * l_prime + 3)/2, y_grid**2/(y_grid**2 + 1))
        
        if n_r > 0 and k < n_r:
            G_2 += spc.gamma(l + gamma + k + 2) * spc.poch(-n_r + 1, k)/spc.poch(2 * gamma + 1, k) * 2**k/spc.factorial(k) * y_grid**(l + 1)/(y_grid**2 + 1)**((l + gamma + k +2)/2) \
            * spc.hyp2f1((l + gamma + k + 2)/2, (l - gamma - k)/2, (2 * l + 3)/2, y_grid**2/(y_grid**2 + 1))
            F_2 += spc.gamma(l_prime + gamma + k + 2) * spc.poch(-n_r + 1, k)/spc.poch(2 * gamma + 1, k) * 2**k/spc.factorial(k) * y_grid**(l_prime + 1)/(y_grid**2 + 1)**((l_prime + gamma + k +2)/2) \
            * spc.hyp2f1((l_prime + gamma + k + 2)/2, (l_prime - gamma - k)/2, (2 * l_prime + 3)/2, y_grid**2/(y_grid**2 + 1))

    G = M_norm_noimag * C_plus * ((N - kappa) * G_1 - n_r * G_2)

    F = R_norm_noimag * C_minus * ((N - kappa) * G_1 + n_r * G_2)

    dist = G**2 + F**2

    return dist
# %%
if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt

    p_grid = np.linspace(0, 20, 100)

    dist = calc_momentum_distribution_function_rel(1, -1, 1, p_grid)

    print(type(dist))

    plt.plot(p_grid, dist)
    # %%
    plt.show()

    # %%
    dist_func = calc_momentum_distribution_function_rel(1, -1, 1, p_grid)
# %%