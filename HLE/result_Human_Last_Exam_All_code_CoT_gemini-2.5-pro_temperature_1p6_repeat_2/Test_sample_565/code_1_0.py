import networkx as nx
from itertools import combinations

def check_graphs(graphs):
    """
    Checks that all graphs in the list are vertex-transitive
    and that no two graphs in the list are isomorphic.
    """
    for g in graphs:
        if not nx.is_vertex_transitive(g):
            return False
            
    for g1, g2 in combinations(graphs, 2):
        if nx.is_isomorphic(g1, g2):
            return False
            
    return True

def calculate_vt_graph_counts():
    """
    Calculates the number of non-isomorphic vertex-transitive graphs
    on 8 vertices for each degree from 0 to 7.
    """
    counts = [0] * 8
    
    # j=0: Null graph N8
    graphs_j0 = [nx.empty_graph(8)]
    if check_graphs(graphs_j0):
        counts[0] = len(graphs_j0)
        
    # j=7: Complete graph K8
    graphs_j7 = [nx.complete_graph(8)]
    if check_graphs(graphs_j7):
        counts[7] = len(graphs_j7)
        
    # j=1: 4 * K2
    graphs_j1 = [nx.disjoint_union_all([nx.path_graph(2) for _ in range(4)])]
    if check_graphs(graphs_j1):
        counts[1] = len(graphs_j1)

    # j=6: Complement of j=1 graph
    graphs_j6 = [nx.complement(g) for g in graphs_j1]
    if check_graphs(graphs_j6):
        counts[6] = len(graphs_j6)

    # j=2: C8 and 2 * C4
    graphs_j2 = [
        nx.cycle_graph(8),
        nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4))
    ]
    if check_graphs(graphs_j2):
        counts[2] = len(graphs_j2)
        
    # j=5: Complements of j=2 graphs
    graphs_j5 = [nx.complement(g) for g in graphs_j2]
    if check_graphs(graphs_j5):
        counts[5] = len(graphs_j5)

    # j=3: 5 cubic vertex-transitive graphs
    # 1. Two disjoint K4
    g3_1 = nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4))
    # 2. Cube graph Q3
    g3_2 = nx.cubical_graph()
    # 3. Circulant graph C8(1,4)
    g3_3 = nx.circulant_graph(8, [1, 4])
    # 4. Cayley graph of Z4 x Z2
    g3_4 = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(2)]
    g3_4.add_nodes_from(nodes)
    for i in range(4):
        for j in range(2):
            v1 = (i, j)
            # Generators: (1,1), (3,1), (0,1)
            v2 = ((i + 1) % 4, (j + 1) % 2)
            v3 = ((i + 3) % 4, (j + 1) % 2)
            v4 = (i, (j + 1) % 2)
            g3_4.add_edge(v1, v2)
            g3_4.add_edge(v1, v3)
            g3_4.add_edge(v1, v4)
    # 5. Cayley graph of D4
    g3_5 = nx.Graph()
    nodes_d4 = [(i, b) for i in range(4) for b in range(2)]
    g3_5.add_nodes_from(nodes_d4)
    # D4 multiplication: (i,b)*(j,c) = (i + (-1)^b*j mod 4, b+c mod 2)
    # Generators s=(0,1), sr=(3,1), r^2=(2,0)
    s1, s2, s3 = (0, 1), (3, 1), (2, 0)
    gen_list = [s1, s2, s3]
    for i in range(4):
        for b in range(2):
            v = (i, b)
            for j_gen, c_gen in gen_list:
                v_neighbor = ((i + ((-1)**b) * j_gen) % 4, (b + c_gen) % 2)
                g3_5.add_edge(v, v_neighbor)
    
    graphs_j3 = [g3_1, g3_2, g3_3, g3_4, g3_5]
    if check_graphs(graphs_j3):
        counts[3] = len(graphs_j3)

    # j=4: Complements of j=3 graphs
    graphs_j4 = [nx.complement(g) for g in graphs_j3]
    if check_graphs(graphs_j4):
        counts[4] = len(graphs_j4)

    # Final result in the required format
    print(f"[{counts[0]}, {counts[1]}, {counts[2]}, {counts[3]}, {counts[4]}, {counts[5]}, {counts[6]}, {counts[7]}]")

calculate_vt_graph_counts()