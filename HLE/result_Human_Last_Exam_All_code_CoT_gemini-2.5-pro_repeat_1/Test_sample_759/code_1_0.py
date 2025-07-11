# The problem is to find the smallest number of edges 'e' for a simple,
# connected graph with an automorphism group of size 3.
# Based on graph theory, the minimum number of vertices for such a graph is 9.
# For a connected graph with 9 vertices, the number of edges must be at least 8.
# Furthermore, because the automorphism group is C3 and the minimal graph has no
# fixed points, the number of edges must be a multiple of 3.
# The smallest multiple of 3 that is >= 8 is 9.
# The following code describes the construction of such a graph with 9 edges.

# Define the vertices and the base edges for the construction
num_vertices = 9
base_edges = [(0, 1), (0, 3), (0, 4)]
all_edges = set()

# The automorphism is i -> (i + 3) mod 9. We generate all edges by applying
# this automorphism to the base edges.
print("The graph is constructed on 9 vertices {0, 1, ..., 8}.")
print("The edges are generated from 3 base edges under the automorphism i -> (i+3) mod 9.")
print("-" * 30)

# Orbit of (0,1)
orbit1 = []
for k in range(3):
    u, v = (0 + 3*k) % num_vertices, (1 + 3*k) % num_vertices
    edge = tuple(sorted((u, v)))
    all_edges.add(edge)
    orbit1.append(edge)
print(f"Base edge (0,1) generates orbit: {orbit1[0]}, {orbit1[1]}, {orbit1[2]}")

# Orbit of (0,3)
orbit2 = []
for k in range(3):
    u, v = (0 + 3*k) % num_vertices, (3 + 3*k) % num_vertices
    edge = tuple(sorted((u, v)))
    all_edges.add(edge)
    orbit2.append(edge)
print(f"Base edge (0,3) generates orbit: {orbit2[0]}, {orbit2[1]}, {orbit2[2]}")

# Orbit of (0,4)
orbit3 = []
for k in range(3):
    u, v = (0 + 3*k) % num_vertices, (4 + 3*k) % num_vertices
    edge = tuple(sorted((u, v)))
    all_edges.add(edge)
    orbit3.append(edge)
print(f"Base edge (0,4) generates orbit: {orbit3[0]}, {orbit3[1]}, {orbit3[2]}")

# The total number of edges is the size of the set of all generated edges.
smallest_e = len(all_edges)

print("-" * 30)
print("The total number of edges in the final graph is the sum of edges in these orbits.")
print(f"{len(orbit1)} + {len(orbit2)} + {len(orbit3)} = {smallest_e}")
print("\nThe smallest number of edges e is:")
print(smallest_e)
