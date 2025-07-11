import collections

# A memoization cache for the recursive Hosoya Z index calculation
hosoya_memo = {}

def get_graph_key(edges):
    """Generates a canonical, hashable key for a graph based on its edges."""
    # Sort each edge tuple, then sort the list of edges to ensure consistency
    return tuple(sorted(tuple(sorted(edge)) for edge in edges))

def calculate_hosoya_z(edges):
    """
    Recursively calculates the Hosoya Z index of a graph defined by a list of edges.
    The Hosoya Z index is the total number of non-adjacent edge sets (matchings).
    """
    # Convert list to tuple to make it hashable for the memoization cache
    edges = tuple(edges)
    
    # Use a canonical key for the memoization cache
    key = get_graph_key(edges)
    if key in hosoya_memo:
        return hosoya_memo[key]

    # Base case: A graph with no edges has one matching (the empty set).
    if not edges:
        return 1

    # Recursive step based on the formula: Z(G) = Z(G-e) + Z(G-{u,v})
    # Pick the first edge e = (u, v)
    e = edges[0]
    u, v = e
    
    # Term 1: Z(G-e) -> Graph with the edge 'e' removed
    edges_minus_e = edges[1:]
    term1 = calculate_hosoya_z(list(edges_minus_e))
    
    # Term 2: Z(G-{u,v}) -> Graph with edge 'e' and its vertices (u,v) removed
    # This means removing all other edges connected to u or v as well.
    edges_after_node_removal = []
    for edge_prime in edges_minus_e:
        if u not in edge_prime and v not in edge_prime:
            edges_after_node_removal.append(edge_prime)
    term2 = calculate_hosoya_z(edges_after_node_removal)
    
    # Store result in cache and return
    result = term1 + term2
    hosoya_memo[key] = result
    return result

def calculate_zagreb_m1(nodes, edges):
    """
    Calculates the Zagreb M1 index.
    It's the sum of the squares of the degrees of the vertices.
    """
    degrees = collections.defaultdict(int)
    for u, v in edges:
        degrees[u] += 1
        degrees[v] += 1
    
    m1_index = sum(d**2 for d in degrees.values())
    return m1_index

def main():
    """
    Main function to define the molecule and compute the indices.
    """
    # The target molecule is Aspartic Acid.
    # We represent its heavy-atom skeleton as a graph.
    # Atoms (vertices) are numbered 0 through 8.
    # Bonds (edges) connect these numbered atoms.
    aspartic_acid_nodes = list(range(9))
    aspartic_acid_edges = [
        (0, 1), (1, 2), (1, 3), # First carboxyl group and link to alpha-carbon
        (3, 4),                 # Alpha-carbon to Nitrogen
        (3, 5),                 # Alpha-carbon to Beta-carbon
        (5, 6),                 # Beta-carbon to Gamma-carbon
        (6, 7), (6, 8)          # Gamma-carbon (second carboxyl) to its oxygens
    ]

    # Calculate the Zagreb(1) index
    zagreb_m1_index = calculate_zagreb_m1(aspartic_acid_nodes, aspartic_acid_edges)

    # Calculate the Hosoya Z index
    hosoya_z_index = calculate_hosoya_z(aspartic_acid_edges)

    # Calculate the final ratio as per the problem description
    final_ratio = (2 * hosoya_z_index) / zagreb_m1_index

    # Print the final result including each number in the equation
    print("Molecule: Aspartic Acid")
    print(f"Zagreb(1) Index (M1): {zagreb_m1_index}")
    print(f"Hosoya Z Index (Z): {hosoya_z_index}")
    print("\nCalculation: (2 * Z) / M1")
    print(f"2 * {hosoya_z_index} / {zagreb_m1_index} = {final_ratio}")

if __name__ == "__main__":
    main()