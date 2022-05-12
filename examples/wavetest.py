#%%
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])

# %%
def A4(tgrid):
    return np.sin(440 * tgrid)

def E5(tgrid):
    return np.sin(659.25 * tgrid)

def C5(tgrid):
    return np.sin(523.25 * tgrid)

def undertone(prime, other):
    return prime + other
# %%
tgrid = np.linspace(0, 1e-1, 10000)

A4_vals = A4(tgrid)

E5_vals = E5(tgrid)

C5_vals = C5(tgrid)

undertone_vals = undertone(A4_vals, C5_vals)
# %%
plt.plot(tgrid, A4_vals, label='A4')

plt.plot(tgrid, C5_vals, label='C5')

plt.plot(tgrid, undertone_vals, label='undertone')

plt.scatter(tgrid[undertone_vals >= np.max(undertone_vals) - 0.02], undertone_vals[undertone_vals >= np.max(undertone_vals) - 0.02], label='cut')

plt.legend()
# %%
tune = 0.015
period = tgrid[undertone_vals >= np.max(undertone_vals) - tune][1] - tgrid[undertone_vals >= np.max(undertone_vals) - tune][0]
# %%
freq = 1/period
# %%
