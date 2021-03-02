import networkx as nx
import re
import numpy as np
import pandas as pd
from osgeo import ogr
import json
import copy

from utils.coordinate import DutchRDtoWGS84, calc_distance

def simplify_graph(G, n):
    if 'time' in G.nodes[n]:
        return G
    # expand temporally
    is_temp_ext = True
    S = []
    T = []
    p = n
    while is_temp_ext:
        # reset flag
        is_temp_ext = False
        # successors of the current node
        c_nodes = G.successors(p)
        for c in c_nodes:
            # if this node is at the same (spatially projected) vertex
            if c[0] == p[0]:
                is_temp_ext = True
                # combine nodes
                # nx.contracted_nodes(G, p, c, self_loops=False)
                T.append(c)
                p = c
            else:
                # the successor is a spatial node
                S.append(c)
    t_lst = []
    for t in T:
        t_time = G.nodes[t].get('time', [])
        t_lst = t_lst + t_time + [t[1]]
        # t_lst.append(t[1])
        # G = nx.contracted_nodes(G, n, t, self_loops=False)
        for tc in G.successors(t):
            G.add_edge(n, tc)
        for tp in G.predecessors(t):
            if tp != n:
                G.add_edge(tp, n)
        G.remove_node(t)
    nx.set_node_attributes(G, {n:{'time': t_lst}})
    # expand spatially
    for s in S:
        if G.has_node(s):
            G = simplify_graph(G, s)
    return G

def simplify_graph_ext(G, n):
    if 'time' in G.nodes[n]:
        return G
    # expand temporally
    is_temp_ext = True
    S = []
    T = []
    p = n
    while is_temp_ext:
        # reset flag
        is_temp_ext = False
        # successors of the current node
        c_nodes = G.successors(p)
        for c in c_nodes:
            # if this node is at the same (spatially projected) vertex
            if (c[0] == p[0]) & (c[2] == p[2]):
                is_temp_ext = True
                # combine nodes
                # nx.contracted_nodes(G, p, c, self_loops=False)
                T.append(c)
                p = c
            else:
                # the successor is a spatial node
                S.append(c)
    t_lst = []
    for t in T:
        t_time = G.nodes[t].get('time', [])
        t_lst = t_lst + t_time + [t[1]]
        # t_lst.append(t[1])
        # G = nx.contracted_nodes(G, n, t, self_loops=False)
        for tc in G.successors(t):
            G.add_edge(n, tc)
        for tp in G.predecessors(t):
            if tp != n:
                G.add_edge(tp, n)
        G.remove_node(t)
    nx.set_node_attributes(G, {n:{'time': t_lst}})
    # expand spatially
    for s in S:
        if G.has_node(s):
            G = simplify_graph(G, s)
    return G

def traverse_node(G, src, successor_route, Gx):
    # is_circular = False
    if 'pattern' in G.nodes[src]:
        return G, [], G.nodes[src]['route'], G.nodes[src]['pattern']
    # check cyclic
    if src in successor_route:
        # is_circular = True
        return G, [src], [], []
    
    successor_route.append(src)
    
    P = []
    S = {}
    T = []

    # expend this node along the temporal axis
    # also, update possible spatial expansions
    is_temp_ext = True
    s0 = src
    while is_temp_ext:
        # reset flag
        is_temp_ext = False
        # successors of the current node
        suc_node = G.successors(s0)
        for s in suc_node:
            # if this node is at the same (spatially projected) vertex
            if s[0] == s0[0]:
                is_temp_ext = True
                T.append(s)
                s0 = s
            else:
                # the successor is a spatial node
                v = s[0]
                if v in S:
                    S[v].append(s)
                else:
                    S[v] = [s]
    # expand spatially
    Sx = list(Gx.successors(src[0]))
    cir_node_lst = []
    if Sx:
        R = []
        P = []
        # cir_ind = []
        for v in Sx:
            ups_node_lst = S.get(v)
            if not ups_node_lst is None:
                for s in ups_node_lst:
                    # traverse node s
                    G, cir_nodes, Rs, Ps = traverse_node(G, s, successor_route, Gx)
                    cir_node_lst += cir_nodes
                    for r, p in list(zip(Rs, Ps)):#list(zip(G.nodes[s]['route'], G.nodes[s]['pattern'])):
                        r_str = ''.join(str(_) for _ in r)
                        is_new = True
                        for i in range(len(R)):
                            r0 = R[i]
                            if not r0:
                                continue
                            r0_str = ''.join(str(_) for _ in r0)
                            if r_str in r0_str:
                                P[i] += p.copy()
                                # PR_P_dict[r0].append(p)
                                is_new = False
                                # cir_ind.append(i)
                            elif r0_str in r_str:
                                R[i] = r.copy()
                                P[i] += p.copy()
                                is_new = False
                                # cir_ind.append(i)
                        if is_new:
                            R.append(r.copy())
                            P.append(p.copy())
                            # cir_ind.append(len(R))
            else:
                R.append([])
                P.append([])
        for r in R:
            r.insert(0, src[0])#[r.insert(0, src[0]) for r in R]
        P = [p+T+[src] for p in P]
        # add this (route/pattern) information to node attribute
        if not cir_node_lst:
            for s in T+[src]:
                nx.set_node_attributes(G, {s:{'route': R, 'pattern': P}})
    else:
        # finish here
        P=[[src]]
        P[0] = P[0] + T
        R = [[src[0]]]
        # print(P)
        for s in T+[src]:
            nx.set_node_attributes(G, {s:{'route': R, 'pattern': P}})
    try:
        # cir_node_lst.remove(src)
        cir_node_lst = list(filter(lambda x: x != src, cir_node_lst))
    except ValueError:
        pass
    del successor_route[-1]
    return G, cir_node_lst, R, P

def traverse_node_ext(G, src, successor_route, Gx):
    # is_circular = False
    if 'pattern' in G.nodes[src]:
        return G, [], G.nodes[src]['route'], G.nodes[src]['pattern']
    # check cyclic
    if src in successor_route:
        # is_circular = True
        return G, [src], [], []
    
    successor_route.append(src)
    
    P = []
    S = {}
    T = []

    # expend this node along the temporal axis
    # also, update possible spatial expansions
    is_temp_ext = True
    s0 = src
    while is_temp_ext:
        # reset flag
        is_temp_ext = False
        # successors of the current node
        suc_node = G.successors(s0)
        for s in suc_node:
            # if this node is at the same (spatially projected) vertex
            if s[0] == s0[0]:
                is_temp_ext = True
                T.append(s)
                s0 = s
            else:
                # the successor is a spatial node
                v = s[0]
                if v in S:
                    S[v].append(s)
                else:
                    S[v] = [s]
#     # if this is the end of congestion (by judging the last element of the node)
#     if src[2] == 0:
        
    # expand spatially
    Sx = list(Gx.successors(src[0]))
    cir_node_lst = []
    if Sx:
        R = []
        P = []
        # cir_ind = []
        for v in Sx:
            ups_node_lst = S.get(v)
            if not ups_node_lst is None:
                for s in ups_node_lst:
                    # traverse node s
                    G, cir_nodes, Rs, Ps = traverse_node_ext(G, s, successor_route, Gx)
                    cir_node_lst += cir_nodes
                    for r, p in list(zip(Rs, Ps)):#list(zip(G.nodes[s]['route'], G.nodes[s]['pattern'])):
                        r_str = ''.join(str(_) for _ in r)
                        is_new = True
                        for i in range(len(R)):
                            r0 = R[i]
                            if not r0:
                                continue
                            r0_str = ''.join(str(_) for _ in r0)
                            if r_str in r0_str:
                                P[i] += p.copy()
                                # PR_P_dict[r0].append(p)
                                is_new = False
                                # cir_ind.append(i)
                            elif r0_str in r_str:
                                R[i] = r.copy()
                                P[i] += p.copy()
                                is_new = False
                                # cir_ind.append(i)
                        if is_new:
                            R.append(r.copy())
                            P.append(p.copy())
                            # cir_ind.append(len(R))
            else:
                R.append([])
                P.append([])
        for r in R:
            r.insert(0, src[0])#[r.insert(0, src[0]) for r in R]
        P = [p+T+[src] for p in P]
        # add this (route/pattern) information to node attribute
        if not cir_node_lst:
            for s in T+[src]:
                nx.set_node_attributes(G, {s:{'route': R, 'pattern': P}})
    else:
        # finish here
        P=[[src]]
        P[0] = P[0] + T
        R = [[src[0]]]
        # print(P)
        for s in T+[src]:
            nx.set_node_attributes(G, {s:{'route': R, 'pattern': P}})
    try:
        # cir_node_lst.remove(src)
        cir_node_lst = list(filter(lambda x: x != src, cir_node_lst))
    except ValueError:
        pass
    del successor_route[-1]
    return G, cir_node_lst, R, P

def get_prime_pattern_2(G, src_lst, do_beautify=False):
    R_final = []
    P_final = []
    for src in src_lst:
        if 'route' not in G.nodes[src]:
            print(src)
        r_lst = G.nodes[src]['route']
        p_lst = G.nodes[src]['pattern']
        for r, p in list(zip(r_lst, p_lst)):
            r_str = ''.join(str(_) for _ in r)
            is_new = True
            for i in range(len(R_final)):
                r0 = R_final[i]
                if not r0:
                    continue
                r0_str = ''.join(str(_) for _ in r0)
                if r_str in r0_str:
                    # P_final[i] += p
                    P_final[i].update(p)
                    is_new = False
                elif r0_str in r_str:
                    R_final[i] = r
                    # P_final[i] += p
                    P_final[i].update(p)
                    is_new = False
            if is_new:
                R_final.append(r)
                P_final.append(set(p))
    if do_beautify:
        for i in range(len(P_final)):
            p = []
            for r in R_final[i]:
                v_t = [v for v in P_final[i] if v[0] == r]
                v_t.sort()
                p += v_t
            P_final[i] = p
    T_final = []
    for p in P_final:
        t = [24*60*2, -1]
        for n in p:
            t[0]=min(t[0], min(G.nodes[n]['time']+[n[1]]))
            t[1]=max(t[1], max(G.nodes[n]['time']+[n[1]]))
        T_final.append(t)
    return {'route': R_final, 'time': T_final, 'pattern': P_final}

def get_prime_pattern(G, src_lst, do_beautify=False):
    r_lst = [r for s in src_lst for r in G.nodes[s]['route']]
    p_lst = [p for s in src_lst for p in G.nodes[s]['pattern']]

    rp_lst = list(zip(r_lst, p_lst))
    rp_lst = sorted(rp_lst, key=lambda rp: len(rp[0]), reverse=True)

    R_final = []
    P_final = []
    for r, p in rp_lst:
        r_str = ''.join(str(_) for _ in r)
        is_new = True
        for i in range(len(R_final)):
            r0 = R_final[i]
            if not r0:
                continue
            r0_str = ''.join(str(_) for _ in r0)
            if r_str in r0_str:
                # P_final[i] += p
                P_final[i].update(p)
                is_new = False
            elif r0_str in r_str:
                R_final[i] = r
                # P_final[i] += p
                P_final[i].update(p)
                is_new = False
        if is_new:
            R_final.append(r)
            P_final.append(set(p))
    T_final = []
    for p in P_final:
        t = [24*60*2, -1]
        for n in p:
            t[0]=min(t[0], min(G.nodes[n]['time']+[n[1]]))
            t[1]=max(t[1], max(G.nodes[n]['time']+[n[1]]))
        T_final.append(t)
    return {'route': R_final, 'time': T_final, 'pattern': P_final}
        
def gen_congpath_info(pattern_info, link_collection, time_ind_max, time_offset=15, simplify=True, overlap_perc=80):
    cong_path = []
    for clust_id, p_info in enumerate(pattern_info):
        clst_cong_path = []
        for sub_ind, (v_lst, t) in enumerate(zip(p_info['route'], p_info['time'])):
            if len(v_lst) == 1:
                continue
            XY = []
            latlon = []
            lnstr = ogr.Geometry(ogr.wkbLineString)
            wvk_id_list = []
            road_num_list = []
            link_distaces = []
            for v0, v1 in zip(v_lst, v_lst[1:]):
                try:
                    e = link_collection.loc[v0, v1]
                except Exception as e:
                    print(e)
                    print(v0, v1)
                    print(p_info['route'].index(v_lst))
                    return
                wvk_id_list.append(e['WVK_ID'])
                road_num_list.append(e['WEGNUMMER'])
                pts = list(zip(e['DutchRD_X'], e['DutchRD_Y']))
                XY += pts
                link_lnstr = ogr.Geometry(ogr.wkbLineString)
                for p in pts:
                    p_latlon = DutchRDtoWGS84(*p)
                    lnstr.AddPoint_2D(*p_latlon)
                    link_lnstr.AddPoint_2D(*p_latlon)
                    latlon.append(p_latlon)
                link_distaces.append(calc_distance(link_lnstr.ExportToWkt()))
            cong_time_ind = [int(max(1, t[0]-time_offset)), int(min(time_ind_max, t[1]+time_offset))]
            # only keep unique road-number
            road_num_set = set()
            road_num_set_add = road_num_set.add
            road_num = [n for n in road_num_list if not (n in road_num_set or road_num_set_add(n))]
            # check if pattern is big enough
            # > 1 km
            # > 15 mins
            d = calc_distance(lnstr.ExportToWkt())
            if d < 1000:
                continue
            if (t[1] - t[0]) < 15:
                continue
            clst_cong_path.append({'XY': np.array(XY),
                            'latlon': np.array(latlon),
                            'geojson': {
                                'type': 'Feature',
                                'geometry': json.loads(lnstr.ExportToJson()),
                                'properties': []
                            }, 
                            'link': {
                                'wvk_id_list': wvk_id_list,
                                'distance': np.array(link_distaces)
                            },
                            'road_num': road_num,
                            'congtimeindext': np.array(cong_time_ind), 
                            'space_extent': d,
                            'congestion_cluster': {
                                'id': clust_id,
                                'sub_ind': sub_ind
                            }
                            })
        if simplify:
            # eliminate patterns that highly overlap with other patterns
            def calc_overlap_distance(p1, p2):
                link1 = p1['link']
                link2 = p2['link']
                # find overlapping id and distance
                id1 = set(link1['wvk_id_list'])
                id2 = set(link2['wvk_id_list'])
                id_common = set.intersection(id1, id2)
                link1_info = dict(zip(link1['wvk_id_list'], link1['distance']))
                link2_info = dict(zip(link2['wvk_id_list'], link2['distance']))
                link_info = {**link1_info, **link2_info}
                overlap_dist = sum(link_info[wvk] for wvk in id_common)
                # find the longer path
                l1 = np.sum(link1['distance'])
                l2 = np.sum(link2['distance'])
                if l1 > l2:
                    longer_path = 1
                    ratio = overlap_dist / l2
                else:
                    longer_path = 2
                    ratio = overlap_dist / l1
                return ratio, longer_path
            
            def try_add_path(path_list, p, min_ratio):
                if path_list == []:
                    path_list.append(p)
                    return path_list
                for pind, inp in enumerate(path_list):
                    ratio, longer_path = calc_overlap_distance(inp, p)
                    if ratio >= min_ratio:
                        if longer_path == 2:
                            path_list[pind] = p
                            return path_list
                        else:
                            return path_list
                path_list.append(p)
                return path_list
            simp_path = []
            for p in clst_cong_path:
                simp_path = try_add_path(simp_path, p, overlap_perc/100.0)
            clst_cong_path = simp_path
        [cong_path.append(c) for c in clst_cong_path]
    return cong_path