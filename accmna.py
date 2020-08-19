import networkx as nx
import sys
import math
import numpy as np
from MEWC import *
from SPICi import *
import time

node_to_cluster_map = {}
ppi = nx.Graph()
ppi_for_weight = nx.Graph()
sim_ppi = nx.Graph()
sim_ppi_for_merge = nx.Graph()
node_prio = []
NODEs_cliques = {}  # max clique for each node
NODE_weights = {}
clusters = {}
simmap = {}
EDGE_weights = {}
EDGE_weights_sim = {}

# for merge
scores = {}     # (icq,ciq[0],ciq[1],internal_conserved)
merged_to_another = {}


class Candidate:
    m1 = 0
    m2 = 0
    m1k = 0
    m2k = 0
    c = 0
    k = 0
    score = 0.0
    icq = 0.0
    ciqw = 0.0
    ciqe = 0
    cluster = []

class C:
    k = 0
    Ts = 0.9
    Td = 0.9
    ICQ = 0.0
    ICQ = 0.0
    cluster_n = 0
    CIQ_edges = 0
    CIQ_weight = 0.0
    internal_conserved = 0


class M:
    cur_n = 0
    mergedn = 0


def readnetwork(datafile):
    f = open(datafile, 'r')
    lines = f.readlines()
    print(f"datafile opened.")
    k = int(lines[0])  # number of networks
    C.k = k
    print(k)
    print("reading networks file......")
    for i in range(1, k + 1):
        netfile = lines[i].strip()
        print(f"opening  {netfile}  ...")
        try:
            with open(netfile, 'r') as f1:
                nn = 0
                ne = 0
                for j in f1:
                    lineArr = j.strip().split()
                    node1 = lineArr[0]
                    node2 = lineArr[1]
                    if ppi.has_node(node1):
                        node1 = lineArr[0]
                    else:
                        ppi.add_node(node1)
                        sim_ppi.add_node(node1)
                        ppi_for_weight.add_node(node1)
                        NODE_weights[node1] = 0.0
                        node_prio.append(node1)
                        nn = nn + 1
                    if ppi.has_node(node2):
                        node2 = lineArr[1]
                    else:
                        ppi.add_node(node2)
                        sim_ppi.add_node(node2)
                        ppi_for_weight.add_node(node2)
                        NODE_weights[node2] = 0.0
                        node_prio.append(node2)
                        nn = nn + 1
                    if not ppi.has_edge(node1, node2):
                        ppi.add_edge(node1, node2)
                        ppi_for_weight.add_edge(node1, node2, weight=1)
                        ne = ne + 1
                        # node_pair = [node1, node2]
                        # EDGE_weights[node_pair] = 1
                # print(ppi.node())
                print(f"There are {nn} nodes and {ne} edge.")
        except FileNotFoundError:
            print(f"can not find {netfile} ppi file!\n")
    best_weight = {}
    node_Best = {}
    best_sim = {}
    no_self_sim_list = []

    print("reading BLAST similiarity file......")
    for i1 in lines[k + 1:]:
        netfile = i1.strip()
        print(f"opening {netfile}. ")
        try:
            with open(netfile, 'r') as f2:
                for j1 in f2:
                    lineArr = j1.strip().split()
                    sim_node1 = lineArr[0]
                    sim_node2 = lineArr[1]
                    w = float(lineArr[2])
                    # print (sim_node1[0:2])
                    if ppi.has_node(sim_node1):
                        if ppi.has_node(sim_node2):
                            if sim_node1 != sim_node2:
                                sim_ppi.add_edge(sim_node1, sim_node2, weight=w)
                        else:
                            ppi.add_node(sim_node2)
                            sim_ppi.add_node(sim_node2)
                            ppi_for_weight.add_node(sim_node2)
                            NODE_weights[sim_node2] = 0.0
                            node_prio.append(sim_node2)
                            nn = nn + 1
                            if sim_node1 != sim_node2:
                                sim_ppi.add_edge(sim_node1, sim_node2, weight=w)

                    else:
                        ppi.add_node(sim_node1)
                        sim_ppi.add_node(sim_node1)
                        ppi_for_weight.add_node(sim_node1)
                        NODE_weights[sim_node1] = 0.0
                        node_prio.append(sim_node1)
                        nn = nn + 1
                        if ppi.has_node(sim_node2):
                            if sim_node1 != sim_node2:
                                sim_ppi.add_edge(sim_node1, sim_node2, weight=w)

                        else:
                            ppi.add_node(sim_node2)
                            sim_ppi.add_node(sim_node2)
                            ppi_for_weight.add_node(sim_node2)
                            NODE_weights[sim_node2] = 0.0
                            node_prio.append(sim_node2)
                            nn = nn + 1
                            if sim_node1 != sim_node2:
                                sim_ppi.add_edge(sim_node1, sim_node2, weight=w)


                    # visit and find each node's max weight of connected node
                    if best_weight.__contains__(sim_node1):
                        if w > best_weight[sim_node1]:
                            best_weight[sim_node1] = w
                    else:
                        best_weight[sim_node1] = w
                    if best_weight.__contains__(sim_node2):
                        if w > best_weight[sim_node2]:
                            best_weight[sim_node2] = w
                    else:
                        best_weight[sim_node2] = w

                    if not best_sim.__contains__(sim_node1):
                        best_sim[sim_node1] = w
                    elif w > best_sim[sim_node1]:
                        best_sim[sim_node1] = w
                    if not best_sim.__contains__(sim_node2):
                        best_sim[sim_node2] = w
                    elif w > best_sim[sim_node2]:
                        best_sim[sim_node2] = w
                    if sim_node1[0:2] != sim_node2[0:2]:
                        no_self_sim_list.append((sim_node1, sim_node2))
                        nodepair = (sim_node1, sim_node2)
                        simmap[nodepair] = w

                    d1 = {}
                    d2 = {}
                    if node_Best.__contains__(sim_node1):
                        d1 = node_Best[sim_node1]
                    if node_Best.__contains__(sim_node2):
                        d2 = node_Best[sim_node2]
                    if not d1.__contains__(sim_node2[0:2]):
                        d1[sim_node2[0:2]] = 0.0
                    if not d2.__contains__(sim_node1[0:2]):
                        d2[sim_node1[0:2]] = 0.0

                    if d1[sim_node2[0:2]] < w:
                        d1[sim_node2[0:2]] = w
                    if d2[sim_node1[0:2]] < w:
                        d2[sim_node1[0:2]] = w

                    node_Best.setdefault(sim_node1, d1)
                    node_Best.setdefault(sim_node2, d2)
                f2.close()
        except FileNotFoundError:
            print(f"can not find {netfile} BLAST similarity file!\n")
    f.close()
    print("Finished reading BLAST similarities.")
    for tt in no_self_sim_list:
        simmap[tt] = simmap[tt] / (pow(best_sim[tt[0]], 0.5) * pow(best_sim[tt[1]], 0.5))
    no_self_sim_list.clear()
    best_sim.clear()
    # print (best_weight)
    print(f"In total, There are {sim_ppi.number_of_nodes()} nodes and {sim_ppi.number_of_edges()} edges.")
    # filter
    print("Filtering BLAST scores.")
    # print(sim_ppi.edges(data='weight'))
    remove_node = []


    for node1, node2, fil_w in sim_ppi.edges(data='weight'):
        # for node1, nbrs in sim_ppi.adjacency():
        #     for node2, e_data in nbrs.items():
        #         fil_w = e_data['weight']
        if fil_w < beta * node_Best[node1][node2[0:2]] or fil_w < beta * node_Best[node2][node1[0:2]]:
            remove_node.append((node1, node2))
    for node_pair in remove_node:
        sim_ppi.remove_edge(node_pair[0], node_pair[1])
    weight_list = []
    for node1, node2, w in sim_ppi.edges(data='weight'):
        weight_list.append(w)
    x = np.array(weight_list)
    max_w = x.max()
    min_w = x.min()
    ranges = max_w - min_w
    print(f"There are {ppi.number_of_nodes()} ppi nodes and {ppi.number_of_edges()} ppi edges.")
    print(
        f"After all, There are {sim_ppi.number_of_nodes()} sim_ppi nodes and {sim_ppi.number_of_edges()} sim_ppi edges.")


def construct_star(node):
    neighbors = list(sim_ppi.neighbors(node))
    total = 0.0
    matchneighbors_c = [node]
    for ppin_c in neighbors:
        if not node_to_cluster_map.__contains__(ppin_c):
            edgeweight_c = sim_ppi.get_edge_data(node, ppin_c)['weight']
            matchneighbors_c.append(ppin_c)
            total = total + edgeweight_c
    star = nx.Graph(sim_ppi.subgraph(matchneighbors_c))
    # print(star.edges())
    return star


def find_clique(g, rc):
    include = []
    g2 = nx.Graph()
    if rc == -1:
        node = list(g.nodes())[0]
        include.append(node)
    elif rc == 0:
        while g2.number_of_edges() == 0 and len(node_prio) > 0:
            g2.clear()
            if len(node_prio) > 0:
                t_node = ''
                temp_max_w = -1.0
                for u in node_prio:
                    if NODE_weights[u] > temp_max_w:
                        temp_max_w = NODE_weights[u]
                        t_node = u
                # t_node = max(NODE_weights.items(), key=lambda x: x[1])
                g2 = construct_star(t_node).copy()
                if g2.number_of_edges() == 0:
                    node_prio.remove(t_node)
                    # del NODE_weights[t_node]
            else:
                break
            if not g2.number_of_edges() == 0:
                node = list(g2.nodes())
                include.append(node[0])
    else:
        excluded = []
        for node in g.nodes():
            if node_to_cluster_map.__contains__(node):
                excluded.append(node)
        for v in excluded:
            g.remove_node(v)
    clique = []
    if g.number_of_edges() > 0 or g2.number_of_edges() > 0:
        b = []
        if rc == 0:
            #b = find_mewc(g2, C.k, include)
            b = find_cluster(g2)
            clique.extend(b)
            return g, clique
        else:
            #b = find_mewc(g, C.k, include)
            b = find_cluster(g)
            clique.extend(b)
            return g, clique
    else:
        return g, clique


def compute_weight(c):
    result = 0.0
    for u in c:
        for v in c:
            if sim_ppi.has_edge(u, v):
                result += sim_ppi.get_edge_data(u, v)['weight']
    return result / 2.0


def give_overall(cluster):
    icq, kk = calculate_icq(cluster)
    int_cons = 0
    ciq, int_cons = calculate_ciq(cluster)
    if (C.CIQ_edges + ciq[1]) > 0.0:
        score = ((1.0 - alpha) * (C.ICQ * C.cluster_n + icq) / (C.cluster_n + 1)) + (alpha * ((C.CIQ_weight + ciq[0]) /
                                                                                              (C.CIQ_edges + ciq[1])))
    else:
        score = ((1.0 - alpha) * ((C.ICQ * C.cluster_n + icq) / (C.cluster_n + 1)))
    score = round(score, 6)
    result = (icq, ciq[0], ciq[1], int_cons)
    return result, score, kk


# for merge
def give_overall_merge(cluster, cn_change, m_cs, update, cluster_id):
    icq, ww = calculate_icq(cluster)
    ciq, int_cons = calculate_ciq_for_merge(cluster, cluster_id)

    icq_e = 0.0
    ciqe_e = 0.0
    ciqw_e = 0.0
    if m_cs[0] > 0:
        icq_e += scores[m_cs[0]][0]
        ciqw_e += scores[m_cs[0]][1]
        ciqe_e += scores[m_cs[0]][2]
    if isinstance(m_cs[1], int):
        if m_cs[1] > 0:
            icq_e += scores[m_cs[1]][0]
            ciqw_e += scores[m_cs[1]][1]
            ciqe_e += scores[m_cs[1]][2]
    if (C.CIQ_edges + ciq[1] - ciqe_e) > 0:
        score = ((1 - alpha) * ((((C.ICQ * M.cur_n) - icq_e) + icq) / (M.cur_n + cn_change))) + \
                (alpha * ((C.CIQ_weight + ciq[0] - ciqw_e) / (C.CIQ_edges + ciq[1] - ciqe_e)))
    else:
        score = (1 - alpha) * ((((C.ICQ * M.cur_n) - icq_e) + icq) / (M.cur_n + cn_change))
    if update:
        C.ICQ = ((((C.ICQ * M.cur_n) - icq_e) + icq) / (M.cur_n + cn_change))
        C.CIQ_weight = (C.CIQ_weight + ciq[0] - ciqw_e)
        C.CIQ_edges = (C.CIQ_edges + ciq[1] - ciqe_e)
    result = (icq, ciq[0], ciq[1], int_cons)
    return result, score, ww


def calculate_icq(cluster):
    result = 0.0
    sn = 0
    kset_icq = []
    for v in cluster:
        net_name = v[0:2]
        if not kset_icq.__contains__(net_name):
            kset_icq.append(net_name)
        cluster2 = cluster.copy()
        for n in cluster2:
            if v[0:2] != n[0:2]:
                sn += 1
                nodepair1 = (n, v)
                nodepair2 = (v, n)
                if simmap.__contains__(nodepair1):
                    result += simmap[nodepair1]
                elif simmap.__contains__(nodepair2):
                    result += simmap[nodepair2]
                # if sim_ppi.has_edge(v, n):
                #     result = result + sim_ppi.get_edge_data(n, v)['weight']
    if sn != 0:
        result /= sn
    else:
        result = 0.0
    kk = len(kset_icq)
    return result, kk


def calculate_ciq(cluster):
    result = (0.0, 0)
    related_clusters = []
    internal_conserved = 0
    for v in cluster:
        for n in cluster:
            if ppi.has_edge(v, n) and v[0:2] == n[0:2]:
                internal_conserved += 1
        for n in list(ppi.neighbors(v)):
            if node_to_cluster_map.__contains__(n) and v[0:2] == n[0:2]:
                if not related_clusters.__contains__(node_to_cluster_map[n]) and not cluster.__contains__(n):
                    related_clusters.append(node_to_cluster_map[n])
    internal_conserved /= 2
    for i in related_clusters:
        result = add_tuples(result, ciq_between_two(cluster, clusters[i]))
    return result, internal_conserved


# for merge
def calculate_ciq_for_merge(cluster, cluster_id):
    result = (0.0, 0)
    related_clusters = []
    internal_conserved = 0
    for v in cluster:
        for n in cluster:
            if ppi.has_edge(v, n) and v[0:2] == n[0:2]:
                internal_conserved += 1
        for n in list(ppi.neighbors(v)):
            if node_to_cluster_map.__contains__(n) and v[0:2] == n[0:2]:
                if not related_clusters.__contains__(node_to_cluster_map[n]) and not cluster.__contains__(n):
                    related_clusters.append(node_to_cluster_map[n])
    internal_conserved /= 2
    for i_m in related_clusters:
        if i_m < cluster_id:
            result = add_tuples(result, ciq_between_two(cluster, clusters[i_m]))
    return result, internal_conserved


def add_tuples(t1, t2):
    result = (t1[0] + t2[0], t1[1] + t2[1])
    return result


def ciq_between_two(c1, c2):
    mutual_k = 0
    n = 0
    kset_t = []
    for nn in c1:
        kname = nn[0:2]
        if not kset_t.__contains__(kname):
            kset_t.append(kname)
    kset2 = []
    for nn in c2:
        kname = nn[0:2]
        if not kset2.__contains__(kname) and kset_t.__contains__(kname):
            kset2.append(kname)
    mutual_k = len(kset2)
    kset_t.clear()
    # calculate edge weights
    for nn1 in c1:
        for nn2 in c2:
            if ppi.has_edge(nn1, nn2) and nn1[0:2] == nn2[0:2]:
                if not kset_t.__contains__(nn1[0:2]):
                    kset_t.append(nn1[0:2])
                n = n + 1
    e = len(kset_t)
    weight = 0.0
    if e < 2:
        weight = 0.0
    elif mutual_k != 0:
        weight = (e * 1.0) / (mutual_k * 1.0)
    else:
        weight = 0.0
    result = (n * weight, n)
    return result


def update(cluster):
    for v in cluster:
        if node_to_cluster_map.__contains__(v):
            return True
    return False


def update_node_stars(cluster):
    for v in cluster:
        sim_ppi.remove_node(v)
        # node = NODE_weights.pop(v)
    # for v in node_prio:
    #     if not node_to_cluster_map.__contains__(v):
    #         if update(NODEs_cliques[v]):
    #             gg = nx.Graph(construct_star(v))
    #             gg, cl = find_clique(gg, -1)
    #             NODEs_cliques[v] = cl.copy()
    #             NODE_weights[v] = compute_weight(cl)


def construct_neighbour_graph(masters):
    star = nx.Graph()

    neighbors = []
    for v in masters:
        for ppin in list(ppi.neighbors(v)):
            if not node_to_cluster_map.__contains__(ppin):  # exclude clustered proteins
                if not star.has_node(ppin):  # exclude already added nodes
                    if v[0:2] == ppin[0:2]:  # only real ppi neighbours
                        star.add_node(ppin)
                        neighbors.append(ppin)
    neighbors2 = neighbors.copy()
    for p2 in neighbors:
        neighbors2.remove(p2)
        for n2 in neighbors2:
            if p2[0:2] != n2[0:2]:
                if sim_ppi.has_edge(p2, n2):
                    edgeweight = sim_ppi.get_edge_data(p2, n2)['weight']
                    star.add_edge(p2, n2, weight=edgeweight)
    return star


def expand(evaluation_star, ev_star_nodes, v_s):
    if not ev_star_nodes.__contains__(v_s):
        evaluation_star.add_node(v_s)
        ev_star_nodes.append(v_s)
        for n in evaluation_star.nodes():
            if v_s != n:
                weight1 = -1.0
                if v_s[0:2] == n[0:2]:
                    weight1 = 0.0
                elif sim_ppi.has_edge(v_s, n):
                    weight1 = sim_ppi.get_edge_data(v_s, n)['weight']
                if weight1 > 0.0:
                    evaluation_star.add_edge(v_s, n, weight=weight1)
    return evaluation_star, ev_star_nodes


def is_ok(cluster, n):
    for c in cluster:
        weight = -1.0
        if n[0:2] == c[0:2]:
            weight = 0.0
        elif sim_ppi.has_edge(n, c):
            weight = sim_ppi.get_edge_data(n, c)['weight']
        if weight < 0.0:
            return False
    return True


def expand_neighbour_graph(clique):
    backbone = []
    real_backbone = []
    kset_e = []
    result = []
    if 1 < len(clique) < C.k:
        # construct backbone graph
        evaluation_star = nx.Graph()
        ev_star_nodes = []

        # Kset & Backbone
        for v in clique:
            evaluation_star, ev_star_nodes = expand(evaluation_star, ev_star_nodes, v)
            backbone.append(v)
            real_backbone.append(v)
            if not kset_e.__contains__(v[0:2]):
                kset_e.append(v[0:2])
        for cv in real_backbone:
            for v in list(sim_ppi.neighbors(cv)):
                if not kset_e.__contains__(v[0:2]):
                    if not node_to_cluster_map.__contains__(v):
                        if is_ok(real_backbone, v):
                            evaluation_star, ev_star_nodes = expand(evaluation_star, ev_star_nodes, v)
        if evaluation_star.number_of_nodes() > len(backbone):
            to_expand = find_mewc(evaluation_star, C.k, backbone)
            for v in to_expand:
                result.append(v)
        else:
            return clique
    else:
        return clique
    return result

def expand_cluster(clique):
    kset = []
    add_node_can = []
    clique_net = nx.Graph(sim_ppi.subgraph(clique))
    for v in clique:
        if not kset.__contains__(v[0:2]):
            kset.append(v[0:2])
    if len(kset) < C.k:
        for cv in clique:
            for v in list(sim_ppi.neighbors(cv)):
                if not kset.__contains__(v[0:2]):
                    if not node_to_cluster_map.__contains__(v):
                        if not add_node_can.__contains__(v):
                            add_node_can.append(v)
        while len(add_node_can) > 0:
            candidate_support = {}
            nEdge = 0
            for u in add_node_can:
                candidate_support[u] = 0.0
                for n in clique:
                    if sim_ppi.has_edge(u, n):
                        candidate_support[u] += sim_ppi.get_edge_data(u, n)['weight']
            add_node = max(candidate_support, key=candidate_support.get)
            for v in clique:
                if sim_ppi.has_edge(v, add_node):
                    nEdge += 1
            nVertex = len(clique)
            density_t = 2.0 * (clique_net.number_of_edges() + nEdge) / ((nVertex + 1) * nVertex)
            increase = nEdge / (density_t * nVertex)
            if increase >= C.Ts and density_t > C.Td:
                clique.append(add_node)
                add_node_can.remove(add_node)
                clique_net = sim_ppi.subgraph(clique)
            else:
                break
    return clique


def compare_candidates(c1, c2):
    if c1[1] > c2[1]:
        return 1
    elif c1[1] < c2[1]:
        return -1
    else:
        return 0


def extract_backbone(backbonefile):
    print("Backbone Extraction Started.")

    lamda = 0.2
    node_degree = {}
    for n in ppi_for_weight.nodes():
        num = ppi_for_weight.degree(n)
        node_degree[n] = num
    node_degree = sorted(node_degree.items(), key=lambda x: x[1])
    nodes_core = nx.core_number(ppi)

    degree = 0
    while node_degree[0][1] <= 10:
        node_info = node_degree.pop(0)
        node = node_info[0]
        degree = node_info[1]
        if degree == 1:
            neighbours = list(ppi_for_weight.neighbors(node))
            for neighbour in neighbours:
                if ppi_for_weight.has_edge(node, neighbour):
                    NODE_weights[neighbour] = NODE_weights[neighbour] + NODE_weights[node] + \
                                              ppi_for_weight.get_edge_data(node, neighbour)['weight']
        if degree > 1:
            neighbours = list(ppi_for_weight.neighbors(node))
            w_neigbours = 0
            for n in neighbours:
                w_neigbours += ppi_for_weight.get_edge_data(node, n)['weight']
            neighbours1 = neighbours.copy()
            for v1 in neighbours:
                neighbours1.remove(v1)
                for v2 in neighbours1:
                    node_core = nodes_core[node]
                    v1_core = nodes_core[v1]
                    v2_core = nodes_core[v2]
                    if v1 != v2 and v1_core > node_core and v2_core > node_core:
                        N_score = (len(neighbours) * (len(neighbours) - 1)) / 2
                        if ppi_for_weight.has_edge(v1, v2):
                            w = ppi_for_weight.get_edge_data(v1, v2)['weight'] + (
                                        NODE_weights[node] + w_neigbours) / N_score
                            ppi_for_weight.add_edge(v1, v2, weight=w)
                        else:
                            w = (NODE_weights[node] + w_neigbours) / N_score
                            ppi_for_weight.add_edge(v1, v2, weight=w)

    print(f"ppi size:{ppi.number_of_nodes()}")
    for ppin in node_prio:
        sum_weight = 0.0
        sum_weight_sim = 0.0
        for node1 in list(ppi_for_weight.neighbors(ppin)):
            sum_weight += ppi_for_weight.get_edge_data(node1, ppin)['weight']
        for node2 in list(sim_ppi.neighbors(ppin)):
            sum_weight_sim += sim_ppi.get_edge_data(node2, ppin)['weight']
        if len(list(sim_ppi.neighbors(ppin))):
            node_core = nodes_core[ppin]
            t_part1 = NODE_weights[ppin] + lamda * sum_weight
            t_part = node_core*(NODE_weights[ppin] + lamda * sum_weight)
            NODE_weights[ppin] = alpha * t_part + (1 - alpha) * sum_weight_sim / len(list(sim_ppi.neighbors(ppin)))
        else:
            NODE_weights[ppin] = alpha * t_part + (1 - alpha) * sum_weight_sim

    # for prio in sim_ppi.nodes():
    #     gg = nx.Graph(construct_star(prio))
    #     # print(gg.number_of_nodes())
    #     gg, cl = find_clique(gg, -1)
    #     # print(gg.number_of_nodes())
    #     NODEs_cliques[prio] = cl.copy()
    #     NODE_weights[prio] = compute_weight(cl)
    #     # print ("attention")
    #     # print (cl)

    outfile = open(backbonefile, "w")

    i = 0
    C.cluster_n = 0

    candidates = []
    NG = {}

    # add first cluster
    ng1 = nx.Graph()
    ng1, cl1 = find_clique(ng1, 0)
    NG[i] = ng1.copy()

    # scores of first
    result_out = []
    in_sc, ww, cluster_k = give_overall(cl1)
    kandscore = (cluster_k, ww)
    first_c = (0, in_sc, cl1, kandscore)
    candidates.append(first_c)

    while len(candidates) > 0:
        candidate = candidates.pop(0)
        cluster = candidate[2]
        related_cluster = candidate[0]
        str_cluster = []

        if len(cluster) > 0:
            i += 1
            C.cluster_n += 1
            for v in cluster:
                node_prio.remove(v)
                str_cluster.append(v)
                result_out.append(v)
                node_to_cluster_map[v] = i
            result_out.sort()
            for r in result_out:
                outfile.write(r + ' ')
            outfile.write('\n')
            result_out.clear()
            update_node_stars(str_cluster)

            print(f"Cluster {i} has {len(cluster)} nodes. related cluster = {related_cluster}.")
            print(f"Network size: {candidate[3][0]}.")
            print(f"nodeprio size: {len(node_prio)}.")
            print(f"sim_ppi size: {sim_ppi.number_of_nodes()} edge size: {sim_ppi.number_of_edges()}.")

            clusters[i] = cluster.copy()

            # update overall score
            score = candidate[1]
            C.ICQ = ((C.ICQ * (C.cluster_n - 1)) + score[0]) / C.cluster_n
            C.CIQ_weight += score[1]
            C.CIQ_edges += score[2]
            C.internal_conserved += score[3]

            # decide neighbour graph
            g = nx.Graph(construct_neighbour_graph(cluster))
            g, cl_n = find_clique(g, i)
            ngh = expand_cluster(cl_n)
            #ngh = expand_neighbour_graph(cl_n)

            NG[i] = g.copy()

            if len(ngh) > 1:

                sc = (0.0, 0.0, 0, 0)
                inkandscore = (0, 0.0)
                neigh = (i, sc, ngh, inkandscore)
                candidates.append(neigh)
            else:
                NG[i].clear()
            # decide replacement graph
            NG[related_cluster], cl_n2 = find_clique(NG[related_cluster], related_cluster)
            cl2 = expand_cluster(cl_n2)
            #cl2 = expand_neighbour_graph(cl_n2)
            if len(cl2) > 1:
                sc = (0.0, 0.0, 0, 0)
                inkandscore = (0, 0.0)
                first2 = (related_cluster, sc, cl2, inkandscore)
                candidates.append(first2)
            else:
                NG[related_cluster].clear()
            c_size = len(candidates)
            for j in range(0, c_size):
                eval_cand = candidates.pop(0)
                rcn = eval_cand[0]
                if update(eval_cand[2]):
                    NG[rcn], new_cluster_cl = find_clique(NG[rcn], rcn)
                    new_cluster = expand_cluster(new_cluster_cl)
                    #new_cluster = expand_neighbour_graph(new_cluster_cl)
                else:
                    new_cluster = eval_cand[2].copy()
                if len(new_cluster) > 1:
                    score, ww, cluster_k = give_overall(new_cluster)
                    inkandscore = (cluster_k, ww)
                    new_cand = (rcn, score, new_cluster, inkandscore)
                    candidates.append(new_cand)
                else:
                    NG[rcn].clear()
            candidates.sort(key=takevalue, reverse=True)
            # sorted(candidate_for_sort.items(), key=lambda kv: (kv[1], kv[0]))
    outfile.close()


def takevalue(elem):
    return elem[3][1]


def give_cans(i, merge_cans):
    ppi_a = clusters[i]
    a = []
    for v in ppi_a:
        a.append(v)
    exists = False
    aset = []
    c_to_c = []
    c_to_n = []

    for v in a:
        if not aset.__contains__(v[0:2]):
            aset.append(v[0:2])
        for n in list(ppi.neighbors(v)):
            if node_to_cluster_map.__contains__(n):
                if node_to_cluster_map[n] > i:
                    if not merged_to_another.__contains__(node_to_cluster_map[n]):
                        c_to_c.append(node_to_cluster_map[n])
            else:
                c_to_n.append(n)
    for n in c_to_n:
        cross = True
        to_add = []
        for v1 in ppi_a:
            to_add.append(v1)
        for v2 in a:
            if v2[0:2] != n[0:2]:
                if not sim_ppi_for_merge.has_edge(v2, n):
                    cross = False
        if cross and not merged_to_another.__contains__(n):
            to_add.append(n)
            kset = []
            for v3 in to_add:
                if not kset.__contains__(v3[0:2]):
                    kset.append(v3[0:2])
            can = Candidate()
            can.m1 = i
            can.m2 = (-1) * int(n[2:])
            can.m1k = len(aset)
            can.m2k = 1
            can.cluster = to_add
            can.k = len(kset)
            merged = (can.m1, can.m2)
            score, sc, kk = give_overall_merge(to_add, 0, merged, False, C.cluster_n + 1)
            can.score = sc
            can.c = 0
            can.icq = score[0]
            can.ciqw = score[1]
            can.ciqe = score[2]
            merge_cans.append(can)
    for j in c_to_c:
        cross = True
        to_add = []
        for v in ppi_a:
            to_add.append(v)
        b = clusters[j]
        bset = []
        for v in a:
            for n in b:
                if v[0:2] != n[0:2]:
                    if not sim_ppi_for_merge.has_edge(v, n):
                        cross = False
        if cross and not merged_to_another.__contains__(j):
            for v in b:
                to_add.append(v)
                if not bset.__contains__(v[0:2]):
                    bset.append(v[0:2])
            kset = []
            for v in to_add:
                if not kset.__contains__(v[0:2]):
                    kset.append(v[0:2])
            can = Candidate()
            can.m1 = i
            can.m2 = j
            can.m1k = len(aset)
            can.m2k = len(bset)
            can.cluster = to_add
            can.k = len(kset)
            merged = (can.m1, can.m2)
            score, sc, kk = give_overall_merge(to_add, -1, merged, False, C.cluster_n + 1)
            can.score = sc
            can.c = -1
            can.icq = score[0]
            can.ciqw = score[1]
            can.ciqe = score[2]
            merge_cans.append(can)
    return merge_cans


def merge_clusters():
    M.mergedn = 0
    merge_cans = []
    for i in range(C.cluster_n, 0, -1):
        merge_cans = give_cans(i, merge_cans)
    done = []
    while len(merge_cans) > 0:
        cand = Candidate()  # merge_cans.pop()#not finish need change
        cand = max(merge_cans, key=lambda x: x.score)
        merge_cans.remove(cand)
        if not done.__contains__(cand.m1) and not done.__contains__(cand.m2):
            M.mergedn += 1
            C.cluster_n += 1
            merged_to_another[cand.m1] = True
            merged_to_another[cand.m2] = True
            done.append(cand.m1)
            done.append(cand.m2)
            cluster_m = cand.cluster
            clusters[C.cluster_n] = cluster_m
            for m in cluster_m:
                node_to_cluster_map[m] = C.cluster_n

            cl_to_update = []
            for m in cluster_m:
                for n in list(ppi.neighbors(m)):
                    if node_to_cluster_map.__contains__(n):
                        cl_to_update.append(node_to_cluster_map[n])
            for i in cl_to_update:
                if i != C.cluster_n:
                    if i != cand.m1 or i != cand.m2:  # can change
                        if (int(cand.m2) < 0 and int(cand.m1) < i) or (
                                int(cand.m2) > 0 and i > min(int(cand.m1), int(cand.m2))):
                            tt = (i, 0)
                            score, sc, kk = give_overall_merge(clusters[i], 0, tt, True, i)
                            # ciqscore = (score[1], score[2])
                            # add = (score[0], ciqscore)
                            scores[i] = score
            print(f"{len(cl_to_update)}  updated cluster scores.")
            merged = (cand.m1, cand.m2)
            score, sc, kk = give_overall_merge(cluster, cand.c, merged, True, C.cluster_n)
            M.cur_n += cand.c
            # ciqscore = (score[1], score[2])
            # add = (score[0], ciqscore)
            scores[C.cluster_n] = score

            j = 0
            size = len(merge_cans)
            for i in range(0, size):
                in_can = merge_cans.pop(0)
                if (cand.m1 != in_can.m1) and (cand.m1 != in_can.m2) and (cand.m2 != in_can.m1) and (
                        cand.m2 != in_can.m2):
                    if cl_to_update.__contains__(in_can.m1) or cl_to_update.__contains__(in_can.m2):
                        merged2 = (in_can.m1, in_can.m2)
                        score, sc, kk, = give_overall_merge(in_can.cluster, in_can.c, merged2, False, C.cluster_n + 1)
                        in_can.score = sc
                        in_can.icq = score[0]
                        in_can.ciqw = score[1]
                        in_can.ciqe = score[2]
                        j += 1
                    else:
                        icq_e = 0.0
                        ciqe_e = 0.0
                        ciqw_e = 0.0
                        if int(in_can.m1) > 0:
                            icq_e += scores[in_can.m1][0]
                            ciqw_e += scores[in_can.m1][1]
                            ciqe_e += scores[in_can.m1][2]
                        if isinstance(in_can.m2, int):
                            if int(in_can.m2) > 0:
                                icq_e += scores[in_can.m2][0]
                                ciqw_e += scores[in_can.m2][1]
                                ciqe_e += scores[in_can.m2][2]
                        in_can.score = ((1 - alpha) * (
                                    (((C.ICQ * M.cur_n) - icq_e) + in_can.icq) / (M.cur_n + in_can.c))) + (
                                               alpha * ((C.CIQ_weight + in_can.ciqw - ciqw_e) / (
                                                   C.CIQ_edges + in_can.ciqe - ciqe_e)))
                    merge_cans.append(in_can)
            give_cans(C.cluster_n, merge_cans)
            print(f"{j} candidates updated.")
            print(f"{len(merge_cans)} candidates.")
            print(f"{M.mergedn} merging completed.")
        else:
            print("wrong!")


def min(a, b):
    if a < b:
        return a
    else:
        return b


if __name__ == '__main__':

    starttime = time.time()

    C.ICQ = 0.0
    C.CIQ_weight = 0.0
    C.CIQ_edges = 0
    C.internal_conserved = 0

    datafile = sys.argv[1]
    global alpha
    global beta
    alpha = float(sys.argv[2])
    beta = float(sys.argv[3])
    backbonefile = sys.argv[4]
    clusterfile = sys.argv[5]
    print(f"alpha is {alpha}")
    print(f"beta is {beta}")
    readnetwork(datafile)
    sim_ppi_for_merge = sim_ppi.copy()
    # PPI network and score reading finished
    extract_backbone(backbonefile)

    C.ICQ = 0.0
    C.CIQ_weight = 0.0
    C.CIQ_edges = 0
    C.internal_conserved = 0
    M.cur_n = 0
    node_to_cluster_map.clear()

    for i in range(1, C.cluster_n + 1):
        cluster = clusters[i].copy()
        for v in cluster:
            node_to_cluster_map[v] = i
        tt = (0, 0)
        score, sc, ww = give_overall_merge(cluster, 1, tt, True, C.cluster_n)
        scores[i] = score
        M.cur_n += 1

    merge_clusters()
    fw = open(clusterfile, 'w')
    for i in range(1, C.cluster_n + 1):
        if not merged_to_another.__contains__(i):
            cluster = clusters[i]
            #result_out = []
            cluster.sort()
            for v in cluster:
                fw.write(v + ' ')
            fw.write('\n')
            # if len(cluster) > 2:
            #     for v in cluster:
            #         result_out.append(v)
            # result_out.sort()
            # if len(result_out) > 0:
            #     for r in result_out:
            #         fw.write(r + ' ')
            #     fw.write('\n')

    endtime = time.time()
    dtime = endtime - starttime
    print(dtime)
