# %%
from atomaton.momentum_space import radial_function
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
# %%
def test_radial_function_rescaled():
    '''
    Check if implemented rescaled radial function yields similar results to hardcoded radial functions for various n and l.
    '''
    Z = 1
    p_grid = np.linspace(0, 5, 1000)

    for n in range(1, 3):
        y_grid = p_grid * n/Z

        for l in range(0, n):
            G_vals_test = radial_function(n, l, Z, y_grid)
            G_vals_ref = radial_reference(n, l, Z, y_grid)

            plt.plot(y_grid, G_vals_test, label='test')
            plt.plot(y_grid, G_vals_ref, label='ref', linestyle='--')
            plt.legend()
            plt.show()

            if np.linalg.norm(G_vals_ref - G_vals_test) < 1e-10:
                print('Test for rescaled radial functions with parameters n = {0} and l = {1} successful!'.format(n, l))
            else:
                print('shiaaat')
                print('norm =', np.linalg.norm(G_vals_ref - G_vals_test))


def radial_reference(n, l, Z, y_grid):
    '''
    Generate hard-coded reference radial functions taken from Bethe, Salpeter.
    '''
    if n == 1 and l == 0:
        G_vals_ref = radial_1_0(Z, n, y_grid)
        print('yo')
    
    elif n == 2 and l == 0:
        G_vals_ref = radial_2_0(Z, n, y_grid)

    elif n == 2 and l == 1:
        G_vals_ref = radial_2_1(Z, n, y_grid)
    
    else: 
        G_vals_ref = np.zeros_like(y_grid)

    return G_vals_ref


def radial_1_0(Z, n, y_grid):
    G_vals_ref = 4 * np.sqrt(2/np.pi) * Z/n * y_grid/((Z/n * y_grid)**2 + 1)**2
    return G_vals_ref

def radial_2_0(Z, n, y_grid):
    G_vals_ref = 32/np.sqrt(np.pi) * Z/n * y_grid * (4 * (Z/n * y_grid)**2 - 1)/(4 * (Z/n * y_grid)**2 + 1)**3
    return G_vals_ref

def radial_2_1(Z, n, y_grid):
    G_vals_ref = 128/np.sqrt(3 * np.pi) * (Z/n * y_grid)**2/(4 * (Z/n * y_grid)**2 + 1)**3
    return G_vals_ref



# %%
if __name__ == '__main__':
    test_radial_function_rescaled()
# %%

# %%
