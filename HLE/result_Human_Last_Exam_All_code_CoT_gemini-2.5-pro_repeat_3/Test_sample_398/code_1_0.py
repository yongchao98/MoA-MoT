import networkx as nx

# This script demonstrates the properties of the graphs discussed in the reasoning.

# Step 1: Define the maximal planar graph G from our counterexample.
# G is a triangular bipyramid on 5 vertices.
# Vertices 0, 1, 2 form a central triangle. Vertex 3 is one apex, vertex 4 is the other.
G = nx.Graph()
G.add_edges_from([(0, 1), (1, 2), (2, 0),  # The central triangle
                  (0, 3), (1, 3), (2, 3),  # Connections to apex 3
                  (0, 4), (1, 4), (2, 4)]) # Connections to apex 4

# The edge 'e' to be added connects the two non-adjacent apexes.
e = (3, 4)

# Step 2: Create the graph G' by adding edge 'e'.
# The resulting graph G' is the complete graph on 5 vertices, K_5.
G_prime = G.copy()
G_prime.add_edge(*e)

# Step 3: Verify the properties of G and G'.
print("--- Analyzing Graph G (Triangular Bipyramid) ---")
n = G.number_of_vertices()
m = G.number_of_edges()
print(f"Number of vertices (n): {n}")
print(f"Number of edges (m): {m}")

# For a maximal planar graph, the number of edges must be 3*n - 6.
maximal_planar_edge_count = 3 * n - 6
print(f"Expected edges for maximal planar graph (3*n - 6): 3 * {n} - 6 = {maximal_planar_edge_count}")
print(f"Is the edge count correct? {m == maximal_planar_edge_count}")
print(f"Is G planar? {nx.is_planar(G)}")

print("\n--- Analyzing Graph G' (G with one more edge) ---")
n_prime = G_prime.number_of_vertices()
m_prime = G_prime.number_of_edges()
print(f"Number of vertices in G': {n_prime}")
print(f"Number of edges in G': {m_prime}")
print(f"Is G' planar? {nx.is_planar(G_prime)}")

print("\n--- Final Conclusion ---")
print("1. G is a maximal planar graph. Adding edge 'e' makes G' non-planar.")
print("2. Since G' can be made planar by removing one edge ('e'), its crossing number is 1.")
print("   This means G' can be drawn with at most one crossing.")
print("3. The drawing is not unique. In our example, to connect vertices 3 and 4,")
print("   the new edge must cross the central triangle (0,1,2).")
print("   We can choose to cross edge (0,1), (1,2), or (2,0).")
print("   Each choice creates a different drawing with one crossing.")
print("Therefore, G' can be drawn with at most one crossing, but this drawing is not unique.")
