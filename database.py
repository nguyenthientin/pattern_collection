import sqlite3
import numpy as np
import json
from datetime import datetime

def db_store_pattern(**kwargs):
    if ('cursor' not in kwargs) and ('db_file' not in kwargs):
        print('specify database file or provide a connection')
        return -1
    do_close = False
    if 'db_file' in kwargs:
        db_file = kwargs.pop('db_file')
        conn = sqlite3.connect(db_file)
        cur = conn.cursor()
        do_close = True
    elif 'cursor' in kwargs:
        cur = kwargs.pop('cursor')
    if 'table' in kwargs:
        table = kwargs.pop('table')
    else:
        print('specify table')
        return -1
    sql_query = 'INSERT INTO ' + table + ' ('
    val_query_fmt = '('
    values = []
    for key,val in kwargs.items():
        sql_query += key + ','
        val_query_fmt += '?,'
        if key in ['speed', 'flow', 'precipitation']:
            if not isinstance(val, list):
                val = val.tolist()
            values.append(json.dumps(val))
        elif key in ['space_resolution', 'time_resolution',
                     'space_extent', 'time_extent', 'total_delay',
                     'no_disturbances']:
            values.append(val)
        elif key in ['linestring', 'time', 'road_num', 'feature_vector']:
            values.append(json.dumps(val))
        elif key in ['date']:
            values.append(val.strftime('%Y-%m-%d'))
    sql_query = sql_query[:-1] + ') VALUES ' + val_query_fmt[:-1] + ')'
    # print(sql_query)
    # Insert the pattern
    cur.execute(sql_query, values)
    if do_close:
        cur.close()
        conn.close()
        return 1