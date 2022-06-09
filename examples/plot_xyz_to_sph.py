#%%
from scipy.special import sph_harm
from scipy.constants import pi
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
figname = 'xyz_to_sph'
figtitle = 'Transformation of a sphere from Cartesian to spherical polar coordinates'
savedir = '../data/' + figname

# %%
import numpy as np
import matplotlib.pyplot as plt
import pytransform3d.coordinates as pc
from pytransform3d.plot_utils import make_3d_axis


thetas, phis = np.meshgrid(np.linspace(0, np.pi, 11),
                           np.linspace(-np.pi, np.pi, 21))
rhos = np.ones_like(thetas)
spherical_grid = np.column_stack((rhos.reshape(-1), thetas.reshape(-1), phis.reshape(-1)))
cartesian_grid = pc.cartesian_from_spherical(spherical_grid)

ax = make_3d_axis(ax_s=1, unit="m", n_ticks=6)
ax.scatter(cartesian_grid[:, 0], cartesian_grid[:, 1], cartesian_grid[:, 2])
ax.plot(cartesian_grid[:, 0], cartesian_grid[:, 1], cartesian_grid[:, 2])
plt.show()
# %%
