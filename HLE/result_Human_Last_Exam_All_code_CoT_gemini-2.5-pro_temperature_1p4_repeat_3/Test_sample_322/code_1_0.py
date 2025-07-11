# The Three Utilities Problem is equivalent to drawing the K3,3 graph.
# Let's define its properties.
V = 6  # Number of vertices (3 houses + 3 utilities)
E = 9  # Number of edges (3 houses * 3 utilities)

print(f"Step 1: Define the graph's properties.")
print(f"Number of Vertices (V): {V}")
print(f"Number of Edges (E): {E}\n")

# Step 2: Assume the graph is planar and apply Euler's Formula (V - E + F = 2)
# to find the number of faces (F) the drawing would have.
print("Step 2: Assume the graph is planar and use Euler's Formula (V - E + F = 2) to calculate the number of faces (F).")
# F = E - V + 2
F = E - V + 2
print(f"Solving for F: F = {E} - {V} + 2 = {F}")
print(f"So, a planar drawing of this graph must have {F} faces.\n")

# Step 3: Use the properties of this specific graph.
# K3,3 is a bipartite graph, so it has no odd-length cycles.
# This means every face must be bounded by at least 4 edges.
min_edges_per_face = 4
print("Step 3: Use the property that for a bipartite graph, each face must be bounded by at least 4 edges.")
print(f"Minimum edges per face: {min_edges_per_face}\n")

# Step 4: Check for a contradiction.
# Let's count the total number of edge boundaries in the graph.
# Each of the E edges borders exactly two faces. So, the total number of boundaries is 2 * E.
total_edge_boundaries = 2 * E

# From the faces' perspective, the total number of boundaries must be at least F * min_edges_per_face.
min_face_boundaries = F * min_edges_per_face

print("Step 4: We will check if the condition 2*E >= 4*F holds true.")
print("The total number of edge boundaries in any planar graph is exactly 2 * E.")
print(f"2 * E = 2 * {E} = {total_edge_boundaries}")

print("\nSince each of the F faces must have at least 4 edges, the total number of edge boundaries must be at least 4 * F.")
print(f"4 * F = 4 * {F} = {min_face_boundaries}")

print("\nTherefore, for the graph to be planar, the following inequality must be true:")
print(f"2 * E >= 4 * F")
print(f"  {total_edge_boundaries}  >=  {min_face_boundaries}")

# Final conclusion
print("\nStep 5: The conclusion.")
if total_edge_boundaries >= min_face_boundaries:
    print("The inequality is true. A solution might exist (this won't be reached for K3,3).")
else:
    print(f"The inequality {total_edge_boundaries} >= {min_face_boundaries} is FALSE.")
    print("This is a logical contradiction. Our initial assumption that the graph can be drawn on a plane without crossing lines must be incorrect.")
    print("Therefore, it is mathematically impossible to solve the Three Utilities Problem.")
