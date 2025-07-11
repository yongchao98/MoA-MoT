import networkx as nx
from itertools import combinations

def solve_subgraph_counting():
    """
    This script solves the problem by constructing the Gosset graph and the specified subgraph (HoG ID 50698),
    and then computationally finds the number of occurrences of the subgraph within the larger graph.
    """

    # Step 1: Construct the target subgraph H (HoG ID 50698).
    # This graph is the Kneser graph K(8,2).
    # Its vertices are the 28 pairs of elements from an 8-element set,
    # and two vertices are adjacent if the pairs are disjoint.
    H = nx.Graph()
    set_for_H = set(range(8))
    nodes_H = list(combinations(set_for_H, 2))
    H.add_nodes_from(nodes_H)

    for u, v in combinations(nodes_H, 2):
        if not set(u).intersection(set(v)):
            H.add_edge(u, v)

    # Step 2: Construct the Gosset graph G.
    # Its 56 vertices can be represented as the 28 pairs (2-subsets) and
    # 28 co-pairs (6-subsets) of an 8-element set.
    G = nx.Graph()
    set_for_G = set(range(8))

    # To distinguish the two types of vertices, we use a prefix ('p' for pair, 's' for 6-subset).
    v2_nodes = [('p', c) for c in combinations(set_for_G, 2)]
    v6_nodes = [('s', c) for c in combinations(set_for_G, 6)]
    G.add_nodes_from(v2_nodes)
    G.add_nodes_from(v6_nodes)

    # Adjacency rules for the Gosset graph:
    # 1. Two 'p' nodes {a,b} and {c,d} are adjacent if disjoint.
    for u, v in combinations(v2_nodes, 2):
        if not set(u[1]).intersection(set(v[1])):
            G.add_edge(u, v)

    # Helper to get the 2-subset complement of a 6-subset node.
    v6_to_comp = {s_node: frozenset(set_for_G.difference(s_node[1])) for s_node in v6_nodes}

    # 2. Two 's' nodes are adjacent if their 6-subsets intersect in 4 elements,
    # which is equivalent to their 2-subset complements being disjoint.
    for u, v in combinations(v6_nodes, 2):
        comp_u = v6_to_comp[u]
        comp_v = v6_to_comp[v]
        if not comp_u.intersection(comp_v):
            G.add_edge(u, v)

    # 3. A 'p' node {a,b} and an 's' node are adjacent if |{a,b} intersect complement_of_s| == 1.
    for u in v2_nodes:
        for v in v6_nodes:
            set_u = set(u[1])
            comp_v = v6_to_comp[v]
            if len(set_u.intersection(comp_v)) == 1:
                G.add_edge(u, v)

    # Step 3: Find and count the subgraphs using NetworkX's isomorphism matcher.
    # This algorithm can be computationally intensive.
    gm = nx.algorithms.isomorphism.GraphMatcher(G, H)

    # A single subgraph can be found via multiple mappings if H has automorphisms.
    # We count the number of unique sets of vertices in G that form an isomorphic subgraph.
    found_subgraphs = set()
    for mapping in gm.subgraph_isomorphisms_iter():
        # The values of the mapping dict are the nodes in G.
        # We use a frozenset of these nodes to represent the unique subgraph.
        subgraph_nodes = frozenset(mapping.values())
        found_subgraphs.add(subgraph_nodes)
    
    count = len(found_subgraphs)

    print(f"The graph with HoG ID 50698 is the Kneser graph K(8,2).")
    print(f"The number of subgraphs isomorphic to K(8,2) contained in the Gosset graph is: {count}")

solve_subgraph_counting()
<<<2>>>