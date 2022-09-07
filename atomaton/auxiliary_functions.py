#%%
import numpy as np
# %%
def rescale_axes(x_grid, y_vals, isovalue=0.99):
    '''
    Rescales array corresponding to x-axis in a plot such that "isovalue" of the area under the curve is shown in the plot and the x-axis is cut off behind that.
    '''
    # find total area under the curve
    area = np.sum(y_vals)

    # loop over all elements starting from the last one and stop when 1 - isovalue is reached as area
    for ii in range(0, len(y_vals)):
        index = len(y_vals) - ii

        area_index = np.sum(y_vals[index:])

        if area_index/area >= 1 - isovalue:
            rescale_index = index
            break
    
    x_grid_rescaled = x_grid[:rescale_index]
    y_vals_rescaled = y_vals[:rescale_index]

    return x_grid_rescaled, y_vals_rescaled, rescale_index
# %%
