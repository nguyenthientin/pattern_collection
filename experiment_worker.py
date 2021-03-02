from ps_network_graph import gen_nw_graph
from ps_network_cut import nc_heuristic_cut
from ps_network_cut import route_2_matlab
from utils.coordinate import DutchRDtoWGS84
from utils.plots import visualise_routes
from ps_cong_network_traversal import simplify_graph, traverse_node, gen_congpath_info, get_prime_pattern

from ps_cong_network_traversal import simplify_graph_ext

import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import json
import os
from osgeo import ogr
from datetime import datetime, timedelta
import requests
import time
import timeit
from glob import glob
from PIL import Image
import pickle
from sys import platform

from utils.coordinate import DutchRDtoWGS84

import warnings
from skimage.segmentation import chan_vese
from skimage.morphology import label as skimage_label
from utils.preprocess import preprocess_speed
from segmentation.congestion_detection import mark_congested_link
from constants import data_server_url

from route_data import retrieve_parse_route_data
from route_processing import construct_route_info
from pattern_data import retrieve_patterns, retrieve_store_patterns

def check_data(data):
    if len(data['raw']['x']) < 2:
        return False
    return True
def traverse_nw(G):
    Gx = nx.DiGraph()
    # Underlying (spatial) network
    for (v,_,_) in G.nodes:
        Gx.add_node(v)
    for (vi, ti, ci),(vj, tj, cj),_ in G.edges:
        if vi != vj:
            Gx.add_edge(vi,vj)
    # 1. find all source nodes
    in_dgr = dict(G.in_degree)
    src_nd = dict(filter(lambda i: i[1]==0, in_dgr.items()))
    # 2. traverse from each of the sources
    for i, src in enumerate(src_nd):
        # print(i, src, end='\r')
        if not G is None:
            G, cir_lst , r_lst , p_lst = traverse_node(G=G, src=src, successor_route=[], Gx=Gx)
    return G, src_nd
def merge_construct_pattern(G, src_nd):
    return get_prime_pattern(G=G, src_lst=src_nd, do_beautify=True)
def check_pattern(data):
    if len(data['raw']['x']) < 2:
        return False
    if np.sum(np.array(data['speed']) < 65) < (10*10): # 5mins * 1km
        return False
    return True

if __name__ == "__main__":
    t_program_start = timeit.default_timer()
    # shape file
    if platform == 'darwin':
        shpfl = '../patternretrieval/nwb/tmp/test_1.shp'
        work_dir = './test/test_1'
        db_file = './test/test_1/congestion.db'
    else:
    #     shpfl = '../patternretrieval/nwb/01-06-2018/Wegvakken/Highways/Wegvakken_Highway_UniDirection.shp'
        shpfl = '/home/tin/Data/nwb/01-01-2019/Wegvakken/Wegvakken_Highway_UniDirection.shp'
        # shpfl = '../patternretrieval/nwb/small/test_1.shp'
        # shpfl = '../patternretrieval/nwb/tmp/test_1.shp'
    #     work_dir = './experiment/Wegvakken_Highway_UniDirection'
        work_dir = '/home/tin/Data/congestion'
        # work_dir = './experiment/test_1'

        db_file = '/home/tin/Data/congestion/congestion_Feb_2019.db'
        # db_file = './experiment/test_1/database/congestion.db'
    db_table = 'congestion_patterns'

    db_dir = os.path.dirname(db_file)
    if not os.path.exists(db_dir):
        os.mkdir(db_dir)
    """ 1. Partition the network
    """
    print('Partion the network')
    G = gen_nw_graph(shpfl, is_refined=True)
    routes,routes_edge = nc_heuristic_cut(G)
    print('Number of routes: {}'.format(len(routes)))
    """ 2. Retrieve data
    """
    # Preprocess routes to prepare some data such as distances of segments/links along routes, linestrings
    print('Prepare route geojsons')
    route_info, link_collection = construct_route_info(routes_edge, num_process=20)
    
    # Download data for these routes
    print('Download & parse route data')
    def daterange(start_date, end_date):
        for n in range(int((end_date - start_date).days)):
            yield start_date + timedelta(days=n)

    start_date = datetime(2019, 2, 13)
    end_date = datetime(2019, 3, 1)
    from_time = '04:00'
    to_time = '23:00'
    date_fmt = '%Y-%m-%d'

    for date in daterange(start_date, end_date):
        t_process_day_start = timeit.default_timer()

        date_str = date.strftime(date_fmt)
        print('===============================')
        print('Processing {}'.format(date_str))
        date_dir = os.path.join(work_dir, date_str)
        if not os.path.exists(date_dir):
            os.mkdir(date_dir)

        link_collection_complete, retriev_time = retrieve_parse_route_data(
                                            route_info, link_collection, 
                                            from_time, to_time, date_str,
                                            num_thread=30, num_process_get=15,
                                            num_process_parse=5)

        # save
        fn = os.path.join(date_dir, 'link_collection_complete.pkl')
        with open(fn, 'wb') as f:
            pickle.dump(link_collection_complete, f)
            
        fn = os.path.join(date_dir, 'route_data_retrieve_time.pkl')
        with open(fn, 'wb') as f:
            pickle.dump(retriev_time, f)

        fn = os.path.join(date_dir, 'link_collection_complete.pkl')
        with open(fn, 'rb') as f:
            link_collection_complete = pickle.load(f)

        """ 3. Network traversal for congestion patterns
        """
        print('Traverse network for congestion patterns')
        link_collection_complete.set_index(['JTE_ID_BEG', 'JTE_ID_END'], inplace=True)

        # Construct the correponding (congestion) network
        print('Construct congestion network/graph G')
        _any = np.array([l for l in link_collection_complete['link_congestion_any']])
        _all = np.array([l for l in link_collection_complete['link_congestion_all']])
        _up = np.array([l for l in link_collection_complete['link_congestion_upstream']])
        _down = np.array([l for l in link_collection_complete['link_congestion_downstream']])
        G = nx.MultiDiGraph()
        link_cong_mask = _any
        # congestion is marked by 1
        if not np.any(link_cong_mask):
            print('No congestion or No data')
            t_process_day_finish = timeit.default_timer()
            print('Process finished in {:.2f} seconds'.format(t_process_day_finish-t_process_day_start))
            continue
        cong_link_ind, cong_time_ind = np.nonzero(link_cong_mask)
        num_discnt = 0
        for il,it in zip(cong_link_ind, cong_time_ind):
            link_info = link_collection_complete.iloc[il]
            v0 = link_info.name[0]#link_info['JTE_ID_BEG']
            v1 = link_info.name[1]#link_info['JTE_ID_END']
            t0 = it
            t1 = it + 1
            if not _all[il, it]:
                if _up[il, it]:
                    G.add_edge((v0, t0, 1), (v1, t0, 0))
                    G.add_edge((v0, t1, 1), (v1, t1, 0))
                    G.add_edge((v0, t0, 1), (v0, t1, 1))
                    G.add_edge((v1, t0, 0), (v1, t1, 0))
                if _down[il, it]:
                    G.add_edge((v0, t0, 2), (v1, t0, 1))
                    G.add_edge((v0, t1, 2), (v1, t1, 1))
                    G.add_edge((v0, t0, 2), (v0, t1, 2))
                    G.add_edge((v1, t0, 1), (v1, t1, 1))
                if (_up[il, it] == 0) and (_down[il, it] == 0):
                    G.add_edge((v0, t0, 2), (v1, t0, 0))
                    G.add_edge((v0, t1, 2), (v1, t1, 0))
                    G.add_edge((v0, t0, 2), (v0, t1, 2))
                    G.add_edge((v1, t0, 0), (v1, t1, 0))
        #     if discontinuities[il, it]:
        #         num_discnt = num_discnt + 1
        #         G.add_edge((v0, t0, 1), (v1, t0, 0))
        #         G.add_edge((v0, t1, 1), (v1, t1, 0))
        #         G.add_edge((v0, t0, 1), (v0, t1, 1))
        #         G.add_edge((v1, t0, 0), (v1, t1, 0))
                
        #         G.add_edge((v0, t0, 0), (v1, t0, 1))
        #         G.add_edge((v0, t1, 0), (v1, t1, 1))
        #         G.add_edge((v0, t0, 0), (v0, t1, 0))
        #         G.add_edge((v1, t0, 1), (v1, t1, 1))
            else:
                G.add_edge((v0, t0, 1), (v1, t0, 1))
                G.add_edge((v0, t1, 1), (v1, t1, 1))
                G.add_edge((v0, t0, 1), (v0, t1, 1))
                G.add_edge((v1, t0, 1), (v1, t1, 1))

        print('no. nodes {}'.format(G.number_of_nodes()))
        print('no. edges {}'.format(G.number_of_edges()))

        # Simplify G  
        # Combine continuous temporal expansion of a node into one single node
        print('Simplify G')
        # 1. find all source nodes
        in_dgr = dict(G.in_degree)
        src_nd = dict(filter(lambda i: i[1]==0, in_dgr.items()))

        #2. traverse from each of the sources
        for src in src_nd:
            G = simplify_graph_ext(G=G, n=src)
        
        print('no. nodes {}'.format(G.number_of_nodes()))
        print('no. edges {}'.format(G.number_of_edges()))

        # Underlying (spatial) network
        # Traverse congestion network for patterns

        # handle connected components separately
        print('Traverse for (single-route) connected components')
        C = nx.weakly_connected_components(G)
        pattern_info = []
        for i, c in enumerate(C):
            _G = G.subgraph(c)
        #     print('Traversing ...')
            _G, src_nd = traverse_nw(_G)
        #     print('Merging ...')
            info = merge_construct_pattern(_G, src_nd)
        #     %time info = traverse_nw_4_pattern_info(_G)
            pattern_info.append(info)

        # Save pattern_info to file
        fn = os.path.join(date_dir, 'pattern_info.pkl')
        with open(fn, 'wb') as f:
            pickle.dump(pattern_info, f)

        fn = os.path.join(date_dir, 'pattern_info.pkl')
        with open(fn, 'rb') as f:
            pattern_info = pickle.load(f)
        # Construct relevant information for congestion patterns  
        # * geojson
        # * time
        print('Construct space/time info for congestion patterns')
        cong_path = gen_congpath_info(pattern_info, 
                                link_collection_complete, 
                                time_ind_max=link_cong_mask.shape[1], 
                                time_offset=15, overlap_perc=75)

        print('Number of cong paths: {}'. format(len(cong_path)))

        # Save cong_path to file
        fn = os.path.join(date_dir, 'congpath.pkl')
        with open(fn, 'wb') as f:
            pickle.dump(cong_path, f)

        fn = os.path.join(date_dir, 'congpath.pkl')
        with open(fn, 'rb') as f:
            cong_path = pickle.load(f)

        """ 4. Download patterns
        """
        print('Downloading patterns')
        time_offset = datetime.strptime(date.strftime('%Y-%m-%d' + ' ' + from_time), '%Y-%m-%d %H:%M')

        pattern_retriev_info = retrieve_store_patterns(cong_path, 
                                time_offset, date,
                                num_thread=100,
                                db_file=db_file, db_table=db_table)
        fn = os.path.join(date_dir, 'pattern_retriev_info.pkl')
        with open(fn, 'wb') as f:
            pickle.dump(pattern_retriev_info, f)
        t_process_day_finish = timeit.default_timer()
        print('Process finished in {:.2f} seconds'.format(t_process_day_finish-t_process_day_start))
    t_program_finish = timeit.default_timer()
    print('Programe finished in {:.2f} seconds'.format(t_program_finish-t_program_start))
