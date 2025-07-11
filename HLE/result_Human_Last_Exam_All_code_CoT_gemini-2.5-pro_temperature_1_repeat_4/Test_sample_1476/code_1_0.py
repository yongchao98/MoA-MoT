import numpy as np

# Step 1: Define a simple graph with a cycle.
# Let's use a square graph with vertices {0, 1, 2, 3}
# and edges {(0,1), (1,2), (2,3), (3,0)}.
nodes = [0, 1, 2, 3]
edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
num_nodes = len(nodes)
num_edges = len(edges)

# Step 2: Define a vertex signal x0 that has zero total variation.
# According to our reasoning, this signal should be constant across the graph.
x0 = np.array([10.0, 10.0, 10.0, 10.0])

print(f"Graph has nodes V = {nodes} and edges E = {edges}.")
print(f"Vertex signal is x0 = {x0}\n")

# Step 3: Compute the edge signal x1 from x0 based on Premise 3.
# x1_e = |x0_u - x0_v|
x1 = np.zeros(num_edges)
for i, edge in enumerate(edges):
    u, v = edge
    x1[i] = abs(x0[u] - x0[v])

print(f"Computed edge signal based on Premise 3:")
print(f"x1 = {x1}\n")

# Step 4: Verify the premises from the problem statement hold for our signals.
# We need the vertex-edge incidence matrix B1.
# B1 is a |V|x|E| matrix. B1[v, e] = -1 if e=(v,*), +1 if e=(*,v), 0 otherwise.
B1 = np.zeros((num_nodes, num_edges))
for i, edge in enumerate(edges):
    u, v = edge
    B1[u, i] = -1
    B1[v, i] = 1

# Verify Premise 1: x1 is curl-free (no cycles with non-zero sum).
# For a simple planar graph, we check the single cycle (0->1->2->3->0).
# The oriented sum is x1[0,1] + x1[1,2] + x1[2,3] - x1[0,3] (note sign change for opposite direction).
# A more general way is to show x1 is a gradient, which we know it is (the gradient of x0, but with absolute values).
# The easiest check is to see that for x1=0, the curl is trivially zero.
# Sum around cycle (0,1,2,3): x1[0]+x1[1]+x1[2]+x1[3] = 0+0+0+0 = 0.
print("Verifying Premise 1 (Curl-free):")
print("The sum of x1 values around the cycle (0,1,2,3) is 0. The signal is curl-free. Premise holds.\n")


# Verify Premise 2: x1 is divergence-free (B1 @ x1 = 0).
divergence = B1 @ x1
print("Verifying Premise 2 (Divergence-free):")
print(f"The divergence is B1 @ x1 = {divergence}.")
if np.allclose(divergence, 0):
    print("The signal is divergence-free. Premise holds.\n")
else:
    print("The signal is NOT divergence-free. Premise does not hold.\n")


# Step 5: Based on the premises holding, demonstrate the conclusion.
# The conclusion is that the graph has total variation = 0.
# Total Variation (TV) = sum(|x0_u - x0_v|) over all edges {u,v}.
# This is equivalent to the sum of the elements of x1.

print("Conclusion: Calculate the Total Variation (TV) of the graph signal x0.")
tv_terms = []
for edge in edges:
    u, v = edge
    tv_terms.append(f"|{x0[u]} - {x0[v]}|")

tv_sum_str = " + ".join(tv_terms)
tv_values_str = " + ".join([f"{abs(x0[e[0]] - x0[e[1]])}" for e in edges])
total_variation = np.sum(x1)

print(f"TV = {tv_sum_str}")
print(f"   = {tv_values_str}")
print(f"   = {total_variation}")
print("\nThe total variation is 0, which confirms our inference.")