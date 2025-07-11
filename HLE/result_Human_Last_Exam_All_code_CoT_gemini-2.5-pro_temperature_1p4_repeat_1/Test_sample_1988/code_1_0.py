import networkx as nx
from itertools import combinations

def construct_target_graph():
    """
    Constructs the graph with HoG ID 50698, which is the Kneser graph K(8, 2).
    Vertices are edges of K8, represented as 2-element tuples.
    Two vertices are adjacent if the edges are disjoint.
    """
    H = nx.Graph()
    nodes = list(combinations(range(8), 2))
    H.add_nodes_from(nodes)
    for u, v in combinations(nodes, 2):
        if set(u).isdisjoint(set(v)):
            H.add_edge(u, v)
    return H

def construct_gosset_graph():
    """
    Constructs a 56-vertex graph often identified as the Gosset graph.
    It consists of two partitions (V1 and V2), each inducing a K(8, 2) graph.
    Vertices in V1 are tuples (i, j).
    Vertices in V2 are strings f"({i},{j})'" to distinguish them.
    """
    G = nx.Graph()
    base_nodes = list(combinations(range(8), 2))

    # Vertices of type 1: tuples e.g., (0, 1)
    v1 = base_nodes
    # Vertices of type 2: strings e.g., "(0,1)'"
    v2 = [f"{c}'" for c in base_nodes]
    
    G.add_nodes_from(v1)
    G.add_nodes_from(v2)

    # Add edges within each partition to form two K(8,2) graphs
    for partition_nodes in [v1, v2]:
        for u, v in combinations(partition_nodes, 2):
            # Parse nodes back to tuples for set operations
            u_tup = eval(u[:-1]) if isinstance(u, str) else u
            v_tup = eval(v[:-1]) if isinstance(v, str) else v
            if set(u_tup).isdisjoint(set(v_tup)):
                G.add_edge(u, v)

    # Add edges between the two partitions
    for u_v1 in v1:
        for v_v2 in v2:
            # v_v2 is a string, parse it back to a tuple
            v_tup = eval(v_v2[:-1])
            # Adjacent if they share exactly one vertex in K8
            if len(set(u_v1).intersection(set(v_tup))) == 1:
                G.add_edge(u_v1, v_v2)
                
    return G

def find_subgraphs():
    """
    Constructs the graphs and finds the number of induced subgraphs.
    Note: The search can be computationally intensive.
    """
    G = construct_gosset_graph()
    H = construct_target_graph()
    
    num_edges_H = H.number_of_edges()

    # The GraphMatcher finds non-induced subgraphs. We will filter them.
    gm = nx.isomorphism.GraphMatcher(G, H)
    
    found_subgraphs_nodesets = set()
    
    # Iterate over all subgraph isomorphisms
    iso_iter = gm.subgraph_isomorphisms_iter()
    for mapping in iso_iter:
        # A mapping's keys form the node set of the subgraph in G.
        nodeset = frozenset(mapping.keys())
        
        # Check if we have already processed this node set
        if nodeset not in found_subgraphs_nodesets:
            # To be an *induced* subgraph, the edge count must match.
            subgraph_induced_by_G = G.subgraph(nodeset)
            if subgraph_induced_by_G.number_of_edges() == num_edges_H:
                found_subgraphs_nodesets.add(nodeset)
    
    print(f"The number of subgraphs with HoG graph ID 50698 contained in the Gosset graph is: {len(found_subgraphs_nodesets)}")

if __name__ == '__main__':
    find_subgraphs()
<<<2>>>