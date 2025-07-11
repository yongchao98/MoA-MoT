import networkx as nx

def analyze_graph_crossing():
    """
    Analyzes the addition of an edge to a maximal planar graph (octahedron)
    to determine if the 1-crossing drawing is unique.
    """
    # A maximal planar graph G is a planar graph where adding any edge makes it non-planar.
    # The resulting graph G' = G U {e} will have a crossing number of 1.
    # We investigate if the drawing of G' is unique by checking how many different edges
    # the new edge 'e' can cross.

    # Step 1: Create an example of a maximal planar graph, the octahedron.
    # n=6 vertices, 3n-6 = 12 edges.
    G = nx.Graph()
    # Let vertex 0 be the "north pole" and 5 be the "south pole".
    # Vertices 1, 2, 3, 4 form the "equator".
    edges = [
        (0, 1), (0, 2), (0, 3), (0, 4),  # North pole to equator
        (5, 1), (5, 2), (5, 3), (5, 4),  # South pole to equator
        (1, 2), (2, 3), (3, 4), (4, 1)   # Edges on the equator
    ]
    G.add_edges_from(edges)

    # Step 2: Choose a non-existent edge 'e'. Let's connect the poles.
    u, v = 0, 5
    e = (u, v)

    print(f"Consider the octahedron graph G, which is maximal planar.")
    print(f"We want to add the edge e = {e} to connect the two non-adjacent 'poles'.")
    print("-" * 40)
    
    # Step 3: To draw e with one crossing, it must cross an existing edge.
    # The candidate edges to cross are those that connect the common neighbors of u and v.
    common_neighbors = sorted(list(nx.common_neighbors(G, u, v)))
    
    print(f"The vertices {u} and {v} have the following common neighbors: {common_neighbors}")

    # Step 4: Identify the specific edges that can be crossed. These are the edges
    # on the "equator" that form a cycle separating u and v.
    candidate_edges_to_cross = []
    for i in range(len(common_neighbors)):
        cn1 = common_neighbors[i]
        cn2 = common_neighbors[(i + 1) % len(common_neighbors)] # Next neighbor in cycle
        if G.has_edge(cn1, cn2):
            candidate_edges_to_cross.append(tuple(sorted((cn1, cn2))))

    print(f"To create a drawing of G U {e} with one crossing, the new edge {e}")
    print("can be drawn to cross any one of the following 'equator' edges:")
    for edge in candidate_edges_to_cross:
        print(f" - Edge {edge}")

    print("-" * 40)
    # Step 5: Conclude based on the number of candidates.
    num_candidates = len(candidate_edges_to_cross)
    print(f"There are {num_candidates} different edges that can be chosen for the crossing.")
    
    if num_candidates > 1:
        print("Since there is more than one choice, the drawing with one crossing is NOT unique.")
    else:
        print("The drawing with one crossing is unique.")
    
    print("\nThis demonstrates that G' can be drawn with at most one crossing, but this drawing is not unique.")


if __name__ == "__main__":
    analyze_graph_crossing()