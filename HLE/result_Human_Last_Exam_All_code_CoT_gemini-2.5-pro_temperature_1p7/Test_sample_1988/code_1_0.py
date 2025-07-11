import networkx as nx
from urllib import request
import itertools

def solve_graph_containment():
    """
    This function calculates the number of subgraphs with HoG ID 50698
    (Kneser graph K(8,2)) contained in the Gosset graph.
    """
    # Step 1: Define the target graph, Kneser graph K(8,2).
    # Its vertices are the C(8,2) edges of a K8, and are adjacent if disjoint.
    G_target = nx.kneser_graph(8, 2)
    num_nodes_G = G_target.number_of_nodes()
    num_edges_G = G_target.number_of_edges()
    print(f"The target subgraph (Kneser graph K(8,2)) has {num_nodes_G} vertices and {num_edges_G} edges.")

    # Step 2: Define the host graph, the Gosset graph.
    # We construct it from the adjacency list on A.E. Brouwer's website.
    url = "https://www.win.tue.nl/~aeb/graphs/sr/Gosset.adj"
    H_gosset = nx.Graph()
    try:
        response = request.urlopen(url)
        # File uses 1-based indexing, so we convert to 0-based.
        # Format is "node: neighbor1 neighbor2 ..."
        adj_data = response.read().decode('utf-8').splitlines()[1:]
        for line in adj_data:
            parts = line.split(':')
            if len(parts) != 2: continue
            
            u = int(parts[0].strip()) - 1
            neighbors_str = parts[1].strip()
            
            if not H_gosset.has_node(u):
                H_gosset.add_node(u)

            if neighbors_str:
                for v_str in neighbors_str.split():
                    v = int(v_str) - 1
                    H_gosset.add_edge(u, v)
    except Exception as e:
        print(f"Error: Could not retrieve or parse the Gosset graph data from {url}.")
        print(f"Details: {e}")
        return

    num_nodes_H = H_gosset.number_of_nodes()
    num_edges_H = H_gosset.number_of_edges()
    print(f"The host graph (Gosset graph) has {num_nodes_H} vertices and {num_edges_H} edges.")

    # Step 3 & 4: Find and count unique induced subgraphs.
    print("\nSearching for induced subgraphs isomorphic to K(8,2)...")
    matcher = nx.isomorphism.GraphMatcher(H_gosset, G_target)
    
    found_subgraphs_nodes = set()
    
    # Iterate through all isomorphism mappings.
    for mapping in matcher.subgraph_isomorphisms_iter():
        # A mapping dictionary's values are the nodes of the subgraph in the host graph.
        subgraph_nodes_in_H = list(mapping.values())
        
        # An isomorphism is for an INDUCED subgraph if the number of edges is identical.
        # The iterator already ensures all edges from the pattern exist in the subgraph.
        # If the edge count is the same, no extra edges can be present.
        if H_gosset.subgraph(subgraph_nodes_in_H).number_of_edges() == num_edges_G:
            # We found an induced subgraph. Add its vertex set to our collection.
            # Using a frozenset makes the set of nodes hashable for adding to a set.
            found_subgraphs_nodes.add(frozenset(subgraph_nodes_in_H))

    # Step 5: Report the final count.
    final_count = len(found_subgraphs_nodes)
    print("\nFinal count:")
    print(f"The number of subgraphs with HoG ID 50698 found in the Gosset graph is {final_count}.")

if __name__ == "__main__":
    solve_graph_containment()