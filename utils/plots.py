import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np

def visualise_routes(routes, shape_fn, fig=None):
    if fig is None:
        fig = plt.figure()
    # Read the (backgound) network
    gdf = gpd.read_file(shape_fn)
    gdf.set_index('WVK_ID', inplace=True)
    # Link the route-info outcome with the geo-dataframe by WVK_ID
    WVK_ID_lst = [[e[2]['WVK_ID'] for e in r] for r in routes]

    # Add an axis
    ax = fig.add_axes([0,0,1,1])
    # plot the background
    # gdf.to_crs('EPSG:28992', inplace=True)
    gdf.plot(ax=ax, color='darkgray', zorder=1, lw=1)

    # prepare colormap for all the routes
    num_colors = len(WVK_ID_lst)
    cmap = plt.get_cmap('Accent')
    colors = [cmap(i) for i in np.linspace(0, 1, num_colors)]

    # now, plot the routes
    for wvk, col in zip(WVK_ID_lst, colors):
        route_df = gdf.loc[wvk]
        plot_route(route_df=route_df, ax=ax, color=col, lw=2)
    # return fig

def plot_route(route_df, ax=None, color='blue', lw=1):
    if ax is None:
        ax = plt.subplot()
    # now, plot the route
    # route_df.to_crs('EPSG:28992', inplace=True)
    route_df.plot(ax=ax, color=color, zorder=2, lw=lw)