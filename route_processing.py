from multiprocessing import Pool
import pandas as pd
import numpy as np
from osgeo import ogr
import json

from utils.coordinate import DutchRDtoWGS84

def _calc_route_info(_route_data):
    link_collection = pd.DataFrame(columns=['WVK_ID','WEGNUMMER',
                                            'JTE_ID_BEG','JTE_ID_END','POS_TV_WOL','BST_CODE','geometry',
                                            'DutchRD_X', 'DutchRD_Y', 'latitude', 'longitude',
                                            'route_id', 'link_offset', 'distance_offset'])
    route_id = _route_data['route_id']
    route_link_collection = _route_data['link_collection']
    XY = []
    latlon = []
    wvk_id = []
    link_begin_ind = []
    link_distance = []
    lnstr = ogr.Geometry(ogr.wkbLineString)
    distance_offset = 0
    for link_ind, link in enumerate(route_link_collection):
        link_info = {'route_id': route_id, 'link_offset': link_ind,
                    'JTE_ID_BEG': link[0], 'JTE_ID_END': link[1], 
                    **link[2]}
        # DutchRD X,Y
        link_geo = ogr.CreateGeometryFromJson(link[2]['geometry'])
        pts = link_geo.GetPoints()
        link_info['DutchRD_X'] = np.array([p[0] for p in pts])
        link_info['DutchRD_Y'] = np.array([p[1] for p in pts])
        link_info['distance_offset'] = distance_offset
        XY += pts
        wvk_id.append(link[2]['WVK_ID'])
        link_begin_ind.append(len(latlon) + 1)
        segment_length = np.zeros(len(pts))
        p_pre = np.array(pts[0])
        for p_ind, p in enumerate(pts):
            p_wgs = DutchRDtoWGS84(*p)
            lnstr.AddPoint_2D(*p_wgs)
            latlon.append(p_wgs)
            p = np.array(p)
            segment_length[p_ind] = np.linalg.norm(p - p_pre)
            p_pre = p
        link_distance.append(np.sum(segment_length))
        link_info['segment_distances'] = segment_length
        link_info['distance_along_route'] = np.cumsum(segment_length) + distance_offset
        link_collection = link_collection.append(link_info, ignore_index=True)
        distance_offset = link_info['distance_along_route'][-1]
    route_info = {'id': route_id,
                    'XY': XY,
                    'latlon': latlon,
                    'link': {'wvk_id': wvk_id, 
                             'distance': link_distance,
                             'beginpointind': link_begin_ind},
                    'geojson': {'type': 'Feature',
                                'geometry': json.loads(lnstr.ExportToJson()),
                                'properties':[]}}
    return route_info, link_collection

def construct_route_info(routes_edge, num_process=4):
    route_data = []
    for id, route in enumerate(routes_edge):
        for l in route:
            if not isinstance(l[2]['geometry'], str):
                l[2]['geometry'] = l[2]['geometry'].ExportToJson()
        route_data.append({'route_id': id, 'link_collection': route})

    with Pool(processes=num_process) as pool:
        res = pool.map(_calc_route_info, (route_data))

    route_info = []
    link_collection = pd.DataFrame(columns=['WVK_ID','WEGNUMMER',
                                            'JTE_ID_BEG','JTE_ID_END','POS_TV_WOL','BST_CODE','geometry',
                                            'DutchRD_X', 'DutchRD_Y', 'latitude', 'longitude',
                                            'route_id', 'link_offset', 'distance_offset'])
    for _data in res:
        _route_info = _data[0]
        coords = [{'lat': c[0], 'lng': c[1]} for c in _route_info['geojson']['geometry']['coordinates']]
        _route_info['geojson']['geometry']['coordinates'] = coords
        route_info.append(_route_info)

        link_collection = link_collection.append(_data[1], ignore_index=True)
    
    # link_collection.drop('geometry', axis='columns', inplace=True)
    route_info = sorted(route_info, key=lambda r: sum(r['link']['distance']), reverse=True)
    return route_info, link_collection