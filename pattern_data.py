from queue import Queue
from threading import Thread, Lock
from multiprocessing import Queue as Process_Queue
import time
import requests
import os
from datetime import datetime, timedelta
import timeit
from PIL import Image
import numpy as np
import json
import re
import sqlite3

from constants import data_server_url
from database import db_store_pattern

def _get_total_delay(stats_str):
    m = re.findall('(\S*) =\s*(\d+\.\d+)', stats_str)
    stats = dict(m)
    total_delay = stats['TotalVehicleLossHours']
    return float(total_delay)

def check_pattern(data):
    if len(data['raw']['x']) < 2:
        return False
    if np.sum(np.array(data['speed']) < 60) < (10*10): # 5mins * 1km
        return False
    return True

def _get_data(pattern_info_queue, time_offset, date, data_dir=None, img_dir=None, db_file=None, db_table='congestion_patterns', lock=None, pattern_retriev_info_queue=None):
    post_data = {
        'service': 'linestring2ndw',
        'start_date': date.strftime('%Y-%m-%d'),
        'space_resolution': 100,
        'time_resolution': 30,
        'docrosssecion': 1,
        'doraw': 1,
        'dostats': 1,
        'douserclass': 0,
        'dotraveltime': 0
    }
    while not pattern_info_queue.empty():
        _data = pattern_info_queue.get()
        i = _data['id']
        p = _data['info']
        # print('Retrieving pattern {} ...'.format(i))
        # print(i, end=' ')
        # pt_start = timeit.default_timer()
        coords = [{'lat': c[0], 'lng': c[1]} for c in p['geojson']['geometry']['coordinates']]
        p['geojson']['geometry']['coordinates'] = coords
        ct_begin = time_offset + timedelta(seconds=int(p['congtimeindext'][0]*60))
        ct_end = time_offset + timedelta(seconds=int(p['congtimeindext'][1]*60))
        post_data['linestring'] = p['geojson']
        post_data['from_time'] = ct_begin.strftime('%H:%M')
        post_data['to_time'] = ct_end.strftime('%H:%M')
        mirror_data_server = 'http://mirrors-ndw.citg.tudelft.nl/service'
        try:
            res = requests.post(data_server_url, json=post_data, verify=False)
        except Exception:
            res = requests.post(mirror_data_server, json=post_data, verify=False)
        try:
            data = res.json()['data']
            if not check_pattern(data):
                pattern_retriev_info_queue.put({'id': i, 'info': 'ignore since too small'})
                pattern_info_queue.task_done()
                continue
            if data_dir:
                pfn = os.path.join(data_dir, '{}.txt'.format(i))
                with open(pfn, 'w') as f:
                    json.dump(data, f)
            if img_dir:
                speed = np.array(data['speed'])
                max_speed = np.max(speed, initial=-1, where=speed<200)
                speed[speed == 99999] = max_speed
                speed = np.uint8(speed/max_speed * 255)
                im = Image.fromarray(speed, 'L')
                im.save(os.path.join(img_dir, '{}.jpg'.format(i)))
            # write to dabase
            if db_file:
                db_data = {}
                db_data['speed'] = data['speed']
                db_data['flow'] = data['flow']
                db_data['space_resolution'] = post_data['space_resolution']
                db_data['time_resolution'] = post_data['time_resolution']
                db_data['linestring'] = post_data['linestring']
                db_data['time'] = [post_data['from_time'], post_data['to_time']]
                db_data['road_num'] = p['road_num']
                db_data['date'] = date
                db_data['space_extent'] = p['space_extent']
                db_data['time_extent'] = int((ct_end - ct_begin).total_seconds()/60)
                db_data['total_delay'] = _get_total_delay(data['stats'])
                # lock.acquire()
                with lock:
                    with sqlite3.connect(db_file) as conn:
                        cursor = conn.cursor()
                        db_store_pattern(cursor=cursor, 
                                        table=db_table,
                                        **db_data)
                        conn.commit()
                # lock.release()
        except Exception as e:
            pattern_retriev_info_queue.put({'id': i, 'info': 'error no data'})
        # pt_finish = timeit.default_timer()
        # print('finished in {:.2f} seconds'.format(pt_finish-pt_start))
        pattern_info_queue.task_done()
    # print('finish thread')

def retrieve_patterns(congestion_paths, time_offset, date, work_dir, num_thread):
    pattern_dir = os.path.join(work_dir, 'patterns')
    if not os.path.exists(pattern_dir):
        os.mkdir(pattern_dir)
    data_dir = os.path.join(pattern_dir, 'text')
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    img_dir = os.path.join(pattern_dir, 'img')
    if not os.path.exists(img_dir):
        os.mkdir(img_dir)

    pattern_info_queue = Queue()
    for pid, cp in enumerate(congestion_paths):
        pattern_info_queue.put({'id': pid, 'info': cp})
    threads = [Thread(target=_get_data, args=(pattern_info_queue, time_offset, date, data_dir, img_dir)) for _ in range(num_thread)]
    [t.start() for t in threads]
    pattern_info_queue.join()

def retrieve_store_patterns(congestion_paths, time_offset, date, work_dir=None, num_thread=100, db_file=None, db_table='congestion_patterns'):
    if work_dir:
        pattern_dir = os.path.join(work_dir, 'patterns')
        if not os.path.exists(pattern_dir):
            os.mkdir(pattern_dir)
        data_dir = os.path.join(pattern_dir, 'text')
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        img_dir = os.path.join(pattern_dir, 'img')
        if not os.path.exists(img_dir):
            os.mkdir(img_dir)
    else:
        pattern_dir = None
        data_dir = None
        img_dir = None

    if db_file:
        if not os.path.exists(db_file):
            with sqlite3.connect(db_file) as conn:
                print('Creating a congestion database')
                cursor = conn.cursor()
                cursor.execute('''CREATE TABLE congestion_patterns
                                (id INTEGER PRIMARY KEY, 
                                speed string, flow string,
                                space_resolution string,
                                time_resolution string, 
                                linestring string,
                                time string,
                                road_num string,
                                date string,
                                space_extent integer,
                                time_extent integer,
                                total_delay integer,
                                no_disturbances integer,
                                feature_vector string,
                                precipitation string,
                                meta_data string)''')
                conn.commit()
            
    else:
        cursor = None
    pattern_retriev_info_queue = Queue() #Process_Queue()
    pattern_info_queue = Queue()
    for pid, cp in enumerate(congestion_paths):
        pattern_info_queue.put({'id': pid, 'info': cp})
    
    lock = Lock()
    threads = [Thread(target=_get_data, args=(pattern_info_queue, time_offset, date, data_dir, img_dir, db_file, db_table, lock, pattern_retriev_info_queue)) for _ in range(num_thread)]
    [t.start() for t in threads]
    pattern_info_queue.join()
    pattern_retriev_info_queue.put('end')
    pattern_retriev_info = []
    _info = pattern_retriev_info_queue.get()
    # while not pattern_retriev_info_queue.empty():
    while _info != 'end':
        pattern_retriev_info.append(_info)
        _info = pattern_retriev_info_queue.get()
        # pattern_retriev_info.append(pattern_retriev_info_queue.get())
    return pattern_retriev_info
    # print('')
