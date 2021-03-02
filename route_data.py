from multiprocessing import Process, Pool, Manager
from multiprocessing import Queue as Process_Queue
from threading import Thread
from queue import Queue
import time
import requests
import numpy as np
import pandas as pd

from constants import data_server_url
from segmentation.congestion_detection import mark_congested_link

import sys
if sys.platform == 'darwin':
    sys.path.append('/Users/tinnguyen/Work/DiTTlab/Data_Servers/Pythonserver-v1.1.0/source/')
else:
    sys.path.append('/home/tin/Dittlab/Python_Data_Server/data-python-server/source/')

from workers.linestring.worker import run as linestring_worker

def check_data(data):
    if len(data['raw']['x']) < 2:
        return False
    return True

def _retrieve_data(route_info, from_time, to_time, from_date):
    t_begin = time.time()
    post_data = {
        'service': 'linestring2ndw',
        'from_time': from_time,
        'to_time': to_time,
        'start_date': from_date,
        'space_resolution': 100,
        'docrosssecion': 1,
        'doraw': 1,
        'dostats': 0,
        'douserclass': 0,
        'dotraveltime': 0
    }
    post_data['linestring'] = route_info['geojson']
    retry = True
    long_waiting = False
    num_retries = 2
    while retry:
        try:
            q = Queue()
            linestring_worker(post_data, q)

            # res = requests.post(data_server_url, json=post_data, verify=False, timeout=1800)
            retry = False
        except requests.exceptions.Timeout:
            print('Wating for too long: {}'.format(route_info['id']))
            num_retries -= 1
            if num_retries == 0:
                long_waiting = True
                retry = False
        except Exception as e:
            print('Got exception when trying to retrieve route data')
            print(e)
    try:
        data = q.get()
        data['route_id'] = route_info['id']
    except Exception as e:
        # print(e)
        # print('nodata: {}'.format(route_info['id']))
        data = None
    finally:
        t_end = time.time()
        if long_waiting:
            time_last = 2000
        else:
            time_last = t_end - t_begin
        proc_time = {'route_id': route_info['id'],
                     'time': time_last, 
                     'is_data': data is not None}
        return data, proc_time

# def _get_data_process(route_info, from_time, to_time, from_date, data_queue, time_queue=None):
def _get_data_process(info):
    route_info = info['route_info']
    from_time = info['from_time']
    to_time = info['to_time']
    from_date = info['from_date']
    data_queue = info['data_queue']
    time_queue = info['time_queue']

    data, proc_time = _retrieve_data(route_info, from_time, to_time, from_date)
    data_queue.put(data)
    # return proc_time
    time_queue.put(proc_time)

def get_route_data(data_queue, time_queue, route_collection, from_time, to_time, from_date, num_thread, num_process_get, num_process_parse):
    # route_queue = Queue()
    # for route in route_collection:
    #     route_queue.put(route)
    # threads = [Thread(target=_get_data, args=(route_queue, from_time, to_time, from_date, data_queue, time_queue)) for _ in range(num_thread)]
    # [t.start() for t in threads]
    # [t.join() for t in threads]
    get_data_params = []
    for r in route_collection:
        get_data_params.append({'route_info': r,
                                'from_time': from_time,
                                'to_time': to_time,
                                'from_date': from_date,
                                'data_queue': data_queue,
                                'time_queue': time_queue})
    with Pool(processes=num_process_get) as pool:
        pool.map(_get_data_process, (get_data_params), chunksize=1)


    # route_queue.join()
    print('finish retrieving route data\n')
    for _ in range(num_process_parse):
        data_queue.put('done')

    # return proc_time

def _parse_data(data_queue, parsed_data_queue, link_collection):
    data = data_queue.get()
    while data != 'done':
        # process data
        if data is not None:
            try:
                # print('parsing data on route {}'.format(data['route_id']))
                if check_data(data):
                    link_route_df = link_collection[
                                        link_collection['route_id']==data['route_id']]
                    lcm,_ = mark_congested_link(
                                            np.array(data['speed']),
                                            link_route_df,
                                            np.array(data['x']),
                                            np.array(data['raw']['x']))
                    parsed_data_queue.put(lcm)
                else:
                    parsed_data_queue.put(None)
            except Exception:
                parsed_data_queue.put(None)
        else:
            parsed_data_queue.put(None)
        # get next piece of data
        data = data_queue.get()
    # print('finish process')
            
def parse_route_data(data_queue, num_route, link_collection, num_process=4):
    parsed_data_queue = Process_Queue()
    processes = [Process(target=_parse_data, args=(data_queue, parsed_data_queue, link_collection)) for _ in range(num_process)]
    [p.start() for p in processes]
    # [p.join() for p in processes]
    link_congestion_mark = pd.DataFrame(columns=['link_congestion_all', 
                                                 'link_congestion_any',
                                                 'link_congestion_upstream', 'link_congestion_downstream'])
    # print('Reading parsed data on Routes: ', end=' ')
    for i in range(num_route):
        # print(i, end=' ')
        print('{}/{}'.format(i+1,num_route), end='\r')
        # print('Receiving parsed data on Routes: {}'.format(i))
        try:
            lcm = parsed_data_queue.get()
            if lcm is not None:
                link_congestion_mark = link_congestion_mark.append(lcm)
        except Exception as exc:
            # terminate all process and return
            print('waiting data for too long. Going to terminate all processes')
            [p.terminate() for p in processes]
            break
        
    print('{}/{}'.format(i+1,num_route))
    link_collection_complete = pd.concat([link_collection, link_congestion_mark], axis='columns')
    link_collection_complete.dropna(subset=['link_congestion_all'], inplace=True)
    print('terminating all parsing processes ...', end='')
    [p.terminate() for p in processes]
    print('Done')
    return link_collection_complete

def retrieve_parse_route_data(route_collection, link_collection, from_time, to_time, from_date, num_thread=20, num_process_get=15, num_process_parse=5):
    # data_queue = Process_Queue()
    m = Manager()
    data_queue = m.Queue()
    time_queue = m.Queue()

    # proc_time = get_route_data(data_queue, time_queue, route_collection, from_time, to_time, from_date, num_thread, num_process_get, num_process_parse)
    get_data_process = Process(target=get_route_data, args=(data_queue, time_queue, route_collection, from_time, to_time, from_date, num_thread, num_process_get, num_process_parse), daemon=False)
    get_data_process.start()

    link_collection_complete = parse_route_data(data_queue, len(route_collection), link_collection, num_process=num_process_parse)
    proc_time = []

    print('terminating data retrieving process ...', end='')
    get_data_process.terminate()
    print('Done')
    # get_data_process.terminate()
#     for _ in range(len(route_collection)):
    while not time_queue.empty():
        time_info = time_queue.get()
        proc_time.append(time_info)

    return link_collection_complete, proc_time
