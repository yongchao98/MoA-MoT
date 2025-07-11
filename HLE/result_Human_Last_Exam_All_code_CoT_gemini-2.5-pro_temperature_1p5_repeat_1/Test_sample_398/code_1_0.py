import networkx as nx

# Step 1: Explain the concept and build a concrete example.
print("Step 1: Understanding Maximal Planar Graphs and constructing an example.")
print("A graph G on n vertices is maximal planar if it is planar and adding any edge between non-adjacent vertices makes it non-planar.")
print("For n >= 3, this means every face in a planar drawing is a triangle, and the number of edges is m = 3n - 6.")
print("\nLet's construct an example: a triangular bipyramid on 5 vertices.")

n = 5
G = nx.Graph()
# Vertices for the triangular bipyramid:
# 0, 1, 2 form the "equator"
# 3 is the "north pole", 4 is the "south pole"
G.add_nodes_from(range(n))

# Add equator edges
G.add_edge(0, 1)
G.add_edge(1, 2)
G.add_edge(2, 0)

# Connect north pole (3) to the equator
G.add_edge(3, 0)
G.add_edge(3, 1)
G.add_edge(3, 2)

# Connect south pole (4) to the equator
G.add_edge(4, 0)
G.add_edge(4, 1)
G.add_edge(4, 2)

m = G.number_of_edges()
print(f"Our graph G has n={n} vertices and m={m} edges.")
print(f"Checking the formula for maximal planarity: 3 * n - 6 = 3 * {n} - 6 = {3 * n - 6}.")
print(f"Since m = {m} equals {3*n-6}, G has the required number of edges.")

# Step 2: Verify properties of G.
print("\nStep 2: Verify the properties of G.")
# For NetworkX versions < 3.0, planarity check returns a tuple.
# For versions >= 3.0, it raises an exception if not planar.
try:
    is_planar, _ = nx.check_planarity(G)
except nx.NetworkXException:
    is_planar = False # Should not happen for this G
print(f"Is G planar? {is_planar}")
print("The only pair of non-adjacent vertices in G is (3, 4), the two poles.")
e = (3, 4)

# Step 3: Form G' = G + e and analyze it.
print("\nStep 3: Form G' by adding the edge e=(3,4) and analyze it.")
G_prime = G.copy()
G_prime.add_edge(e[0], e[1])

# Check what G' is
K5 = nx.complete_graph(5)
is_K5 = nx.is_isomorphic(G_prime, K5)
print(f"The resulting graph G' is isomorphic to K5 (the complete graph on 5 vertices): {is_K5}.")

try:
    is_prime_planar, _ = nx.check_planarity(G_prime)
except nx.NetworkXException as ex: # Catches non-planarity for newer networkx
    print(f"Is G' planar? False ({ex})")
    is_prime_planar = False
if is_prime_planar: # For older networkx
     print(f"Is G' planar? {is_prime_planar}")

print("As expected, G' is non-planar because G was maximal planar.")

# Step 4: Explain why the drawing of G' has 1 crossing and is not unique.
print("\nStep 4: Conclusion about the drawing of G'.")
print("We have established that G' is non-planar. A known theorem states that a graph formed by adding one edge to a planar graph has a crossing number of 1.")
print("This means G' can be drawn with exactly one crossing.")
print("\nIs this drawing unique?")
print("Consider a planar drawing of our example G. The 'equator' triangle (0,1,2) separates the 'north pole' (3) from the 'south pole' (4).")
print("To add the edge e=(3,4), we must cross one of the edges of this separating triangle.")
print("We have a choice:")
print("  - We can draw the edge (3,4) to cross edge (0,1). This is one valid drawing.")
print("  - We can draw the edge (3,4) to cross edge (1,2). This is a second, different drawing.")
print("  - We can draw the edge (3,4) to cross edge (2,0). This is a third, different drawing.")
print("Since there are multiple, distinct ways to draw G' with a single crossing, the drawing is not unique.")
print("\nTherefore, the correct statement is: G' can be drawn in the plane with at most one crossing, but this is not unique.")