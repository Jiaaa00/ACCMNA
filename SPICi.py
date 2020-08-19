import networkx as nx
import pandas as pd

class A:
    Ts = 0.5
    Td = 0.5

node_degree = {}
serach_node = []

def find_cluster(g0):
    g = nx.Graph(g0.copy())
    result = []
    serach_node.clear()
    node_degree.clear()
    for n in g.nodes():
        serach_node.append(n)
        neighbours = g.neighbors(n)
        degree = 0.0
        for v in neighbours:
            degree += g.get_edge_data(n, v)['weight']
        node_degree[n] = degree

    if len(g.nodes()) >= 2:
        seed_node = max(node_degree, key=node_degree.get)
        serach_node.remove(seed_node)
        neighborlist = list(g.neighbors(seed_node))
        if neighborlist:
            result = expand(g, seed_node, serach_node)
        else:
            result.append(seed_node)
    return result


def expand(g, seed_node, serach_node):
    cluster = nx.Graph()
    candidate = []
    node_clustered = [seed_node]
    totalWeight = 0
    second_seed = ''
    seed_neighborlist = list(g.neighbors(seed_node))
    candidate.extend(seed_neighborlist)
    max_weight = 0.0
    for i in seed_neighborlist:
        w = g.get_edge_data(seed_node, i)['weight']
        if w > max_weight:
            max_weight = w
    group0 = []
    group1 = []
    group2 = []
    group3 = []
    group4 = []
    for i in seed_neighborlist:
        w = g.get_edge_data(seed_node, i)['weight']
        r = w / max_weight
        if 0.8 < r <= 1:
            group0.append(i)
        elif 0.6 < r <= 0.8:
            group1.append(i)
        elif 0.4 < r <= 0.6:
            group2.append(i)
        elif 0.2 < r <= 0.4:
            group3.append(i)
        elif 0 < r <= 0.2:
            group4.append(i)
    max_degree = 0.0
    if group0:
        for i in group0:
            if node_degree[i] > max_degree:
                max_degree = node_degree[i]
                second_seed = i
    elif group1:
        for i in group1:
            if node_degree[i] > max_degree:
                max_degree = node_degree[i]
                second_seed = i
    elif group2:
        for i in group2:
            if node_degree[i] > max_degree:
                max_degree = node_degree[i]
                second_seed = i
    elif group3:
        for i in group3:
            if node_degree[i] > max_degree:
                max_degree = node_degree[i]
                second_seed = i
    elif group4:
        for i in group4:
            if node_degree[i] > max_degree:
                max_degree = node_degree[i]
                second_seed = i
    serach_node.remove(second_seed)
    candidate.remove(second_seed)
    node_clustered.append(second_seed)
    cluster = g.subgraph(node_clustered)
    for node in g.neighbors(second_seed):
        if node not in node_clustered and not candidate.__contains__(node):
            candidate.append(node)
    totalWeight += g.get_edge_data(seed_node, second_seed)['weight']
    nVertex = 2
    while len(candidate) > 0:
        candidate_support = {}
        nEdge = 0
        for u in candidate:
            candidate_support[u] = 0.0
            for v in node_clustered:
                if g.has_edge(u, v):
                    candidate_support[u] += g.get_edge_data(u, v)['weight']
        t_node = max(candidate_support, key=candidate_support.get)
        for v in node_clustered:
            if g.has_edge(v, t_node):
                nEdge += 1
        density = 2.0 * cluster.number_of_edges() / (nVertex * (nVertex - 1))
        density_t = 2.0 * (cluster.number_of_edges() + nEdge) / ((nVertex+1) * nVertex)
        increase = nEdge / (density_t * nVertex)
        if increase >= A.Ts and density_t > A.Td:
            node_clustered.append(t_node)
            serach_node.remove(t_node)
            candidate.remove(t_node)
            totalWeight += candidate_support[t_node]
            nVertex += 1
            cluster = g.subgraph(node_clustered)
            for node in g.neighbors(t_node):
                if node not in node_clustered and not candidate.__contains__(node):
                    candidate.append(node)
        else:
            break

    return node_clustered

if __name__ == '__main__':
    net = nx.Graph()
    net.add_edges_from((1))