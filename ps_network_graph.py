import networkx as nx
from osgeo import ogr
import geopandas as gpd
from shapely.geometry import LineString
import os.path

def gen_nw_graph(shpfl, is_refined=False):
    """ Generate a graph that represents the given shapefile.
    """
    G = nx.MultiDiGraph()
    drv = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = drv.Open(shpfl, 0) # 0 means read-only, 1 means writeable
    # Check null
    if dataSource is not None:
        layer = dataSource.GetLayer()
        featureCount = layer.GetFeatureCount()
        print('Number of features in {}: {}'.format(shpfl, featureCount))
    else:
        print ('Could not open {}'.format(shpfl))
    # Create an equivalent graph for this shape file
    G = nx.MultiDiGraph()
    num_bidi = 0
    for feature in layer:
        hecto_ltr = feature.GetField('HECTO_LTTR')
        if hecto_ltr in ['a', 'b', 'c', 'd']:
            continue
        RIJRICHTNG = feature.GetField('RIJRICHTNG') # relative direction (jte on driving)
        WVK_ID = feature.GetField('WVK_ID')         # link/segment id
        JTE_ID_BEG = feature.GetField('JTE_ID_BEG') # node id
        JTE_ID_END = feature.GetField('JTE_ID_END') # node id
        WEGNUMMER = feature.GetField('WEGNUMMER')   # road number
        POS_TV_WOL = feature.GetField('POS_TV_WOL') # driving direction
        BST_CODE = feature.GetField('BST_CODE')     # road function/type
        #linestring = feature.GetGeometryRef().GetPoints() # XY coordinate
        geometry = feature.geometry().Clone()       # geometry
        edge_lst = []
        geometry_lst = []
        if not is_refined:
            if RIJRICHTNG is 'H':
                edge_lst.append([JTE_ID_BEG, JTE_ID_END])
                geometry_lst = [geometry]
                # G.add_edge(JTE_ID_BEG, JTE_ID_END,WEGNUMMER=WEGNUMMER)
            elif RIJRICHTNG is 'T':
                edge_lst.append([JTE_ID_END, JTE_ID_BEG])
                # G.add_edge(JTE_ID_END, JTE_ID_BEG,WEGNUMMER=WEGNUMMER)
                # now, reverse the geometry of this link
                lnstr = ogr.Geometry(ogr.wkbLineString)
                pts = geometry.GetPoints()
                for p in pts[::-1]:
                    lnstr.AddPoint_2D(*p)
                geometry_lst = [lnstr]
            else:
                num_bidi += 1
                edge_lst.append([JTE_ID_BEG, JTE_ID_END])
                edge_lst.append([JTE_ID_END, JTE_ID_BEG])
                # G.add_edge(JTE_ID_BEG, JTE_ID_END,WEGNUMMER=WEGNUMMER)
                # G.add_edge(JTE_ID_BEG, JTE_ID_END,WEGNUMMER=WEGNUMMER)
                lnstr = ogr.Geometry(ogr.wkbLineString)
                pts = geometry.GetPoints()
                for p in pts[::-1]:
                    lnstr.AddPoint_2D(*p)
                geometry_lst = [geometry, lnstr]
                # geometry_lst.append(lnstr)
        else:
            edge_lst.append([JTE_ID_BEG, JTE_ID_END])
            geometry_lst = [geometry]

        for e, geo in zip(edge_lst, geometry_lst):
            G.add_edge(e[0], e[1], key=WVK_ID,
                                WEGNUMMER=WEGNUMMER,
                                WVK_ID=WVK_ID,
                                POS_TV_WOL=POS_TV_WOL,
                                BST_CODE=BST_CODE,
                                geometry=geo)
    print('Number of edges: ',G.number_of_edges())
    print('Number of nodes: ',G.number_of_nodes())
    print('Number of bidirectional nodes', num_bidi)
    # return
    return G

def preprocess_geoshape(shpfl_in, shpfl_out=None):
    # read the input shape file
    gdf = gpd.read_file(shpfl_in)
    gdf.set_index('WVK_ID', inplace=True)
    gdf = gdf[~gdf.index.duplicated(keep='first')]
    # refine the geojson according to the following rules
    # 1. if RIJRICHTNG is 'T', swap nodes
    # 2. if RIJRICHTNG is empty, add a reverse link
    WVK_ID_max = gdf.index.max()
    for i in range(gdf.shape[0]):
        link = gdf.iloc[i]
    #     if link['HECTO_LTTR'] in ['a', 'b', 'c', 'd']:
    #         continue
        if link['RIJRICHTNG'] is 'H':
            continue
        elif link['RIJRICHTNG'] is 'T':
            # reverse this link
            coords = list(link['geometry'].coords)
            gdf['geometry'].iloc[i] = LineString(coords[::-1])        
            gdf.loc[link.name, 'JTE_ID_BEG'], gdf.loc[link.name, 'JTE_ID_END'] = \
                gdf.loc[link.name, 'JTE_ID_END'], gdf.loc[link.name, 'JTE_ID_BEG']
        else:
            WVK_ID_max += 1
            gdf.loc[WVK_ID_max] = link
            coords = list(link['geometry'].coords)
            gdf.loc[WVK_ID_max,'geometry'] = LineString(coords[::-1])
            gdf.loc[WVK_ID_max, 'JTE_ID_BEG'], gdf.loc[WVK_ID_max, 'JTE_ID_END'] = \
                gdf.loc[WVK_ID_max, 'JTE_ID_END'], gdf.loc[WVK_ID_max, 'JTE_ID_BEG']
    # save
    if shpfl_out is None:
        hw_nwb_dir = os.path.split(shpfl_in)[0]
        shpfl_out = os.path.join(hw_nwb_dir, 'Wegvakken_Highway_UniDirection.shp')
    gdf.to_file(driver = 'ESRI Shapefile', filename= shpfl_out)

def extract_highway(shpfl_in, shpfl_out=None):
    nl_gdf = gpd.read_file(shpfl_in)
    # For COSI, we only need highway links
    # Item 'WEGBEHSRT': R (National), P (Provincial)
    # Item 'ROUTELTR': A, N
    isNational = (nl_gdf['WEGBEHSRT'] == 'R')
    isProvincial = (nl_gdf['WEGBEHSRT'] == 'P')
    isARoad = (nl_gdf['ROUTELTR'] == 'A')
    isNRoad = (nl_gdf['ROUTELTR'] == 'N')
    hw_gdf = (nl_gdf[isNational | isARoad])

    if not shpfl_out:
        shp_dir = os.path.dirname(shpfl_in)
        shpfl_out = os.path.join(shp_dir, 'Wegvakken_Highway.shp')
    else:
        shp_dir = os.path.dirname(shpfl_out)
    # tmp_dir = os.path.join(NWB_SHP_DIR, 'tmp')
    if not os.path.exists(shp_dir):
        os.makedirs(shp_dir)
    hw_gdf.to_file(driver = 'ESRI Shapefile', filename= shpfl_out)

if __name__ == "__main__":
    # nwb_shpfl = '/Users/tinnguyen/Work/patternretrieval/nwb/01-01-2019/Wegvakken/Wegvakken.shp'
    # hw_shpfl = '/Users/tinnguyen/Work/patternretrieval/nwb/01-01-2019/Wegvakken/Wegvakken_Highway.shp'
    # shpfl_unidirect = '/Users/tinnguyen/Work/patternretrieval/nwb/01-01-2019/Wegvakken/Wegvakken_Highway_UniDirection.shp'
    nwb_shpfl = '/home/tin/Data/nwb/01-01-2019/Wegvakken/Wegvakken.shp'
    hw_shpfl = '/home/tin/Data/nwb/01-01-2019/Wegvakken/Wegvakken_Highway.shp'
    shpfl_unidirect = '/home/tin/Data/nwb/01-01-2019/Wegvakken/Wegvakken_Highway_UniDirection.shp'
    extract_highway(nwb_shpfl, hw_shpfl)
    preprocess_geoshape(shpfl_in=hw_shpfl, shpfl_out=shpfl_unidirect)
    G = gen_nw_graph(shpfl_unidirect, is_refined=True)
    G.edges(data=True)
