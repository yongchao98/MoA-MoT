import networkx as nx
from itertools import combinations

# 1. Define the maximal planar graph (octahedron)
# Vertices are 0 to 5. Let 0 be 'top', 5 be 'bottom', and 1,2,3,4 be the 'equator'.
G = nx.Graph()
G.add_edges_from([
    (0, 1), (0, 2), (0, 3), (0, 4),  # Top to equator
    (5, 1), (5, 2), (5, 3), (5, 4),  # Bottom to equator
    (1, 2), (2, 3), (3, 4), (4, 1)   # Equator edges in a cycle
])

# 2. Define the new edge 'e' between two non-adjacent vertices
u, v = 0, 5
e = (u, v)

# 3. Analyze the properties
print("Let G be the octahedron graph (which is maximal planar).")
print(f"Let e be the edge {e} connecting the 'top' and 'bottom' vertices, which are not adjacent in G.")

# By definition, G' = G + e is non-planar. A theorem states its crossing number is 1.
# G' can be drawn with exactly one crossing.
print("\nThe resulting graph G' = G + e is non-planar and can be drawn with exactly one crossing.")

# 4. Check for uniqueness
# A one-crossing drawing is created by making the new edge 'e' cross an existing edge in G.
# A property of maximal planar graphs ensures that this is always possible by crossing an
# edge that connects two common neighbors of u and v.
# If there is more than one such edge, the drawing is not unique.

# Find common neighbors of u and v
common_neighbors = sorted(list(nx.common_neighbors(G, u, v)))
print(f"The common neighbors of the vertices {u} and {v} are: {common_neighbors}")

# Find which edges connecting these common neighbors can be "swapped" with 'e'
# to produce a planar graph. Each such "swap" or "flip" corresponds to a valid
# way of drawing G' with one crossing.
possible_crossings = []
for w1, w2 in combinations(common_neighbors, 2):
    if G.has_edge(w1, w2):
        # Temporarily create a graph where (w1, w2) is swapped with (u, v)
        G_flipped = G.copy()
        G_flipped.remove_edge(w1, w2)
        G_flipped.add_edge(u, v)
        # If the result is planar, (w1,w2) is a candidate for the crossing
        is_flipped_planar, _ = nx.check_planarity(G_flipped)
        if is_flipped_planar:
            possible_crossings.append((w1, w2))

print("\nThe number of choices for the single crossing corresponds to the number of edges connecting these common neighbors.")
print("The possible edges in G that the new edge e can cross are:")
for cross_edge in possible_crossings:
    print(f"- Edge {cross_edge}")

print(f"\nSince there are {len(possible_crossings)} different edges that the new edge {e} can cross, the one-crossing drawing of G' is not unique.")
print("\nConclusion: The correct statement is that G' can be drawn with at most one crossing, but this drawing is not unique.")

<<<B>>>