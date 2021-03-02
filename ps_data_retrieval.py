import os
from osgeo import ogr
from utils.coordinate import DutchRDtoWGS84
from datetime import datetime
import requests
import json

import sys
sys.path.append('/Users/tinnguyen/Work/DiTTlab/Data_Servers/Pythonserver-v1.0.0/source/workers')

def retrieve_data(**kwargs):
    """ Retrieve data for all given routes
    """
    # default params
    param_default = {
        'data_dir': os.getcwd()
    }
    # params
    data_dir = kwargs.get('data_dir', default=param_default['data_dir'])
    route_lst = kwargs.get('route_info')
    # prepare routes
    for route in route_lst:
        # get linestring
        ln_str = route_2_geojson(route)


def route_2_geojson(route_edge):
    lnstr = ogr.Geometry(ogr.wkbLineString)
    for e in route_edge:
        pts = e[2]['geometry'].GetPoints()
        for p in pts:
            p = DutchRDtoWGS84(*p)
            lnstr.AddPoint_2D(p[1], p[0])
            # lnstr.AddPoint_2D(*p)
    return lnstr.ExportToJson()

def retrieve_ndw_data(ln_str, time):
    # get detectors
    detectors = get_detectors(ln_str, time[0])
    # request
    dtformat = '%Y-%m-%d %H:%M' # yyyy-mm-dd HH:MM
    ndw_req = {'service': 'ndw',
               'fromlane': 1,
               'tolane': 1,
               'agg': 1,
               'space_resolution': 100,
               'time_resolution':100,
               'detectorInfo': detectors['crossSections'],
               'geojson': {'date': time[0].strftime('%Y-%m-%d'),
                           'type': 'LineString', 
                           'geometry': json.loads(ln_str)},
               'detectors': detectors,
               'start': time[0].strftime(dtformat),
               'end': time[1].strftime(dtformat)}
    url = 'http://mirrors-ndw.citg.tudelft.nl:8081/service'
    r = requests.post(url=url, json=ndw_req)
    return r.json()

def get_detectors(route_linestring, time_from):
    route_json = json.loads(route_linestring)
    url = 'http://dittlab-apps.tudelft.nl/apps/webservice/find/detectors'
    dtformat = '%Y-%m-%dT%H' # yyyy-mm-ddTHH:MM
    dtstr = time_from.strftime(dtformat)
    r = requests.post(url = url, json={"route": route_json['coordinates'],
                                       "dateAndTime": dtstr})
    return r.json()
