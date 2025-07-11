import networkx as nx
from itertools import product

def count_slices(G):
    """
    Counts the number of slices N(G) for a graph G.
    It does this by finding K(G), the number of valid 2-edge-colorings,
    and then N(G) = K(G) / 2.
    A coloring is valid if at no vertex are all incident edges the same color.
    """
    num_edges = G.number_of_edges()
    # Use sorted tuples for edges to handle undirected nature consistently
    edges = sorted([tuple(sorted(e)) for e in G.edges()])
    
    valid_colorings_count = 0
    # Iterate through all 2^|E| edge colorings (0 or 1)
    for colors in product([0, 1], repeat=num_edges):
        edge_to_color_map = dict(zip(edges, colors))
        is_valid_coloring = True
        for node in G.nodes():
            # Get colors of edges incident to the current node
            incident_edge_colors = []
            for neighbor in G.neighbors(node):
                edge = tuple(sorted((node, neighbor)))
                incident_edge_colors.append(edge_to_color_map[edge])
            
            # Check if all colors are the same
            if len(set(incident_edge_colors)) == 1:
                is_valid_coloring = False
                break
        
        if is_valid_coloring:
            valid_colorings_count += 1
            
    # N(G) = K(G) / 2
    # The number of valid colorings must be even, N(G) is an integer.
    return valid_colorings_count // 2

def find_m(n, max_v=12):
    """
    Finds M(n), the smallest number of vertices m for which a cubic graph G
    has N(G) as a multiple of n.
    """
    print(f"Searching for M({n})...")
    # List of (name, generator) for cubic graphs, sorted by size
    # Using a pre-defined list of small cubic graphs
    cubic_graphs_to_check = [
        ("K_4", nx.complete_graph(4)),
        ("Prism(6)", nx.circular_ladder_graph(3)),
        ("K_3,3", nx.complete_bipartite_graph(3, 3)),
        ("Cube(8)", nx.cubical_graph())
        # Add more 8-vertex graphs if needed, but Cube(8) is sufficient for this problem.
    ]

    for name, G in cubic_graphs_to_check:
        num_v = G.number_of_nodes()
        print(f"Checking graph {name} with {num_v} vertices.")
        
        n_g = count_slices(G)
        print(f"Calculated N({name}) = {n_g}.")

        if n_g % n == 0:
            print(f"{n_g} is a multiple of {n}. Found M({n}) = {num_v}.\n")
            return num_v
        else:
            print(f"{n_g} is not a multiple of {n}.\n")
            
    return "none" # Not found within the checked graphs

def solve():
    """
    Solves for M(0), M(3), M(5) and prints the result.
    """
    # M(0): Based on theorem that N(G) is always odd for cubic G. Odd numbers are non-zero.
    m0 = "none"
    print("Determining M(0):")
    print("A 'slice' count N(G) for a cubic graph is always an odd integer.")
    print("A multiple of 0 must be 0. Since 0 is even, N(G) can never be 0.")
    print("Therefore, no such graph exists.\n")
    
    # M(3)
    m3 = find_m(3)
    
    # M(5)
    m5 = find_m(5)

    final_answer = f"{m0},{m3},{m5}"
    print(f"The final result for M(0),M(3),M(5) is: {final_answer}")
    return final_answer

if __name__ == '__main__':
    solve()
    # The final answer format for the platform
    # print(f"<<<{solve()}>>>")