import networkx as nx
import urllib.request
from collections import defaultdict

def calculate_modularity():
    """
    Calculates the Newman-Girvan modularity for the American College Football network.
    
    The graph data is fetched from the canonical source to match the provided image.
    The communities are the NCAA football conferences, which correspond to the node
    colors in the image visualization.
    """
    print("### Calculating Modularity for the College Football Network ###\n")

    # The image and problem refer to the canonical American College Football network dataset.
    # The provided edge list is from a different, smaller graph.
    # To match the image, we fetch the standard dataset from Mark Newman's website.
    url = "http://www-personal.umich.edu/~mejn/netdata/football.gml"
    
    try:
        response = urllib.request.urlopen(url)
        gml_data = [line.decode('utf-8') for line in response.readlines()]
    except Exception as e:
        print(f"Failed to download the dataset. Error: {e}")
        print("Please check your internet connection or the URL.")
        return

    # Load the graph from the GML data
    # The 'label='id'' argument ensures node labels are integers from the GML file
    G = nx.parse_gml(gml_data, label='id')

    # Get the total number of edges (m)
    m = G.number_of_edges()
    two_m = 2 * m
    
    print(f"The network has {G.number_of_nodes()} nodes and {m} edges.\n")

    # Group nodes by their conference (community), which is stored in the 'value' attribute
    communities = defaultdict(list)
    for node, data in G.nodes(data=True):
        community_id = data['value']
        communities[community_id].append(node)
    
    # Sort communities by ID for consistent output
    sorted_communities = dict(sorted(communities.items()))

    print("--- Modularity Calculation ---")
    print(f"Formula: Q = \u03A3 [ (l_c / m) - (d_c / 2m)\u00B2 ]")
    print(f"Total edges m = {m}")
    print(f"Total degree sum 2m = {two_m}\n")
    
    total_q = 0.0
    
    # Iterate over each community to calculate its contribution to the modularity
    for c_id, nodes in sorted_communities.items():
        # Create a subgraph containing only the nodes in the current community
        subgraph = G.subgraph(nodes)
        
        # l_c: number of edges within the community
        l_c = subgraph.number_of_edges()
        
        # d_c: sum of degrees of nodes in the community (using the original graph)
        d_c = sum(deg for _, deg in G.degree(nodes))
        
        term1 = l_c / m
        term2 = (d_c / two_m) ** 2
        community_q = term1 - term2
        total_q += community_q
        
        print(f"--- Community {c_id} (Nodes: {len(nodes)}) ---")
        print(f"  Edges within community (l_c): {l_c}")
        print(f"  Sum of degrees (d_c): {d_c}")
        print(f"  Term 1 (l_c / m)         = {l_c} / {m} = {term1:.5f}")
        print(f"  Term 2 (d_c / 2m)\u00B2        = ({d_c} / {two_m})\u00B2 = {term2:.5f}")
        print(f"  Contribution to Q        = {term1:.5f} - {term2:.5f} = {community_q:.5f}\n")

    print("--------------------------------")
    print(f"Final Modularity Q = {total_q:.5f}")
    
    # Verification using the built-in networkx function
    community_list_for_nx = [set(nodes) for nodes in sorted_communities.values()]
    nx_q = nx.algorithms.community.modularity(G, community_list_for_nx)
    print(f"Verification (networkx.modularity) = {nx_q:.5f}")

if __name__ == "__main__":
    calculate_modularity()