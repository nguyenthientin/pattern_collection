import networkx as nx
from ps_network_graph import gen_nw_graph
import json
from osgeo import ogr
from utils.coordinate import DutchRDtoWGS84

def nc_heuristic_cut(G, min_len=10):
    """ Cut a graph, which represents a (highway) network, into non-overlapping routes
    The algorithm is essentially grows routes from multiple sourses.
    Some HEURISTIC rules can be defined:
    1. First extract routes on the same highway (by numbers)
    2. Prioritize main roads when chosing next outgoing links

    Parameters:
    min_len (int): the minimum length of routes when applying the first heuristic rule

    Returns:
    route_lst (list): list of routes
    """
    route_lst = []
    route_edge_lst = []
    # Heuristic 1 - Unique highways
    # Get all highways
    highway_lst = list(set([data['WEGNUMMER'] for _,_,data in G.edges(data=True)]))
    highway_lst = sorted(highway_lst, key=lambda n: 'ZZZ' if n is None else n)
    # directions
    hwdir = ['L', 'R']
    # process each highway
    for hw in highway_lst:
        for d in hwdir:
            # get the sub_graph that represents this highway (hw) on this direction (d)
            e = [(u,v,k) for u,v,k,a in G.edges(data=True, keys=True) 
                                    if a['WEGNUMMER'] == hw and a['POS_TV_WOL'] == d]
            _g = G.edge_subgraph(e).copy()
            # let's cut this sub_network
            rs, rs_edge = nc_randomwalk(_g)
            # keep long routes
            rs = [r for r in rs if len(r) >= min_len]
            rs_edge = [re for re in rs_edge if len(re) >= (min_len-1)] # n_edges = n_nodes - 1
            # check new routes
            if rs:
                # update the final routes
                route_lst += rs
                route_edge_lst += rs_edge
                # remove related edges from the graph
                for re in rs_edge:
                    edges = [e[0:2] + (e[2]['WVK_ID'],) for e in re]
                    G.remove_edges_from(edges)
    # Process the remainder
    # cut the network
    rs, rs_edge = nc_randomwalk(G)
    route_lst += rs
    route_edge_lst += rs_edge
    return route_lst, route_edge_lst

def nc_randomwalk(G):
    """ Cut the graph G into RANDOM non-overlapping routes.
    """
    _g = G.copy()
    is_new_route = True #indicate more routes are possible
    route_lst = []
    route_edge_lst = []
    while is_new_route:
        # find original nodes in_degree = 0
        ori_node_id_lst = [n for n in _g.nodes() if _g.in_degree[n] == 0 and _g.out_degree[n]>0]
        ori_edge_lst = [e for n_id in ori_node_id_lst for e in _g.out_edges(n_id,data=True)]
        if ori_edge_lst:
            # yes, there are new route(s)
            is_new_route = True
        else:
            # no more route (or the remaining is circular)
            # TODO: handle circular
            is_new_route = False
            pass
        # prioritize main road 'HR'
        ori_edge_lst.sort(key=lambda e: e[2]['BST_CODE']=='HR', reverse=True)
        # find/extract all routes starting from these origin nodes
        for e in ori_edge_lst:
            is_expandable = True
            r = [e[0]]
            re = []
            while is_expandable:
                # extend routes
                r.append(e[1])
                # extend (edge-based) routes
                re.append(e)
                # remove the edge
                _g.remove_edge(*e[0:2])
                # check outgoing edges from e[1]
                out_edges = list(_g.out_edges(e[1],data=True))
                if out_edges:
                    # prioritize main road 'HR'
                    out_edges.sort(key=lambda e: e[2]['BST_CODE']=='HR', reverse=True)
                    # pick the first edge
                    e = out_edges[0]
                    # expandable
                    is_expandable = True
                else:
                    is_expandable = False
            # save this route
            route_lst.append(r)
            route_edge_lst.append(re)
    return route_lst, route_edge_lst

def route_2_matlab(route_edge_lst):
    RouteInfo = []
    for re in route_edge_lst:
        XY = []
        latlon = []
        wvk_id = []
        link_begin_ind = []
        lnstr = ogr.Geometry(ogr.wkbLineString)
        for e in re:
            pts = e[2]['geometry'].GetPoints()
            XY += pts
            wvk_id.append(e[2]['WVK_ID'])
            link_begin_ind.append(len(latlon) + 1)
            for p in pts:
                p = DutchRDtoWGS84(*p)
                lnstr.AddPoint_2D(*p)
                latlon.append(p)
        RouteInfo.append({'XY': XY,
                        'latlon': latlon,
                        'link': {'wvk_id': wvk_id, 'beginpointind': link_begin_ind},
                        'geojson': {'type': 'Feature',
                                    'geometry': json.loads(lnstr.ExportToJson()),
                                    'properties':[]}})
    return RouteInfo

if __name__ == "__main__":
    # shpfl = '/Users/tinnguyen/Work/patternretrieval/nwb/01-06-2018/Wegvakken/Highways/Wegvakken_Highway.shp'
    shpfl = '/Users/tinnguyen/Work/patternretrieval/nwb/small/SouthHolland.shp'
    G = gen_nw_graph(shpfl)
    routes,routes_edge = nc_heuristic_cut(G)
    RouteInfo = route_2_matlab(routes_edge)
    with open('RouteInfo.json', 'w') as f:
        json.dump(RouteInfo, f)
    print(len(routes))
