import numpy as np

def print_matrix(name, m):
    """Helper function to print matrices and vectors nicely."""
    print(f"--- {name} ---")
    # check if it's a vector and print it in a more readable way
    if len(m.shape) == 1 or m.shape[1] == 1:
        print(np.round(m.flatten(), 2))
    else:
        print(np.round(m, 2))
    print()

# Step 1: Define the graph G (a 4-cycle)
# V = {0, 1, 2, 3}
# E = {e0=(0,1), e1=(1,2), e2=(2,3), e3=(3,0)}
num_vertices = 4
num_edges = 4

# The graph has no triangles, so B2 is an empty matrix.
# B2 maps triangles to edges. Shape: |E| x |T| = 4x0
B2 = np.zeros((num_edges, 0))

# Define the vertex-edge incidence matrix B1
# We choose an orientation for each edge: 0->1, 1->2, 2->3, 3->0
# B1 has shape |V| x |E|. B1[v, e] is -1 if e starts at v, +1 if e ends at v.
B1 = np.array([[-1,  0,  0,  1],  # Vertex 0 connected to e0(out) and e3(in)
               [ 1, -1,  0,  0],  # Vertex 1 connected to e0(in) and e1(out)
               [ 0,  1, -1,  0],  # Vertex 2 connected to e1(in) and e2(out)
               [ 0,  0,  1, -1]]) # Vertex 3 connected to e2(in) and e3(out)

print_matrix("B1 (Vertex-Edge Incidence Matrix)", B1)

# Step 2: Define a non-trivial signal x0 on vertices
# Let's create a signal that is not constant.
x0 = np.array([0, 10, 0, 10])
print_matrix("x0 (Signal on Vertices)", x0)

# Step 3: Compute the edge signal x1 from x0
# x1_e = |x0_u - x0_v|. This is the magnitude of the gradient of x0.
# The gradient is B1.T @ x0
grad_x0 = B1.T @ x0
x1 = np.abs(grad_x0)
print_matrix("Gradient of x0 (B1.T @ x0)", grad_x0)
print_matrix("x1 (Signal on Edges, |grad_x0|)", x1)

# Step 4: Verify the given conditions for our example
# Condition 1: "no cycles with non-zero sum" -> B2.T @ x1 = 0
# For our graph, B2 is empty, so this is trivially true.
# B2.T has shape 0x4.
curl_x1 = B2.T @ x1
print("Condition 1 check (curl-free): B2.T @ x1")
print(f"Result is a vector of shape {curl_x1.shape}, which is a zero vector. Condition holds.\n")


# Condition 2: "B1 @ x1 = 0"
div_x1 = B1 @ x1
print_matrix("Condition 2 check (divergence-free): B1 @ x1", div_x1)
# We check if the result is close to a zero vector.
if np.allclose(div_x1, 0):
    print("Condition B1 @ x1 = 0 holds.\n")
else:
    print("Condition B1 @ x1 = 0 does not hold.\n")

# Step 5: Check the options
# Option C: x1 is in ker(B1). This is the same as Condition 2, which holds.
# Option D: Total Variation = 0. TV = sum(x1).
tv = np.sum(x1)
print(f"Option D check: Total Variation is sum(x1) = {tv}.")
if tv == 0:
    print("Option D holds.\n")
else:
    print("Option D is false for this non-trivial signal.\n")


# Option A: x1 is an eigenvector of the edge Laplacian L1.
# The prompt defines L1 = B0 @ B0.T + B1.T @ B1.
# B0 is the boundary operator on vertices, so B0=0. Thus L1 = B1.T @ B1.
# (Note: For a graph with triangles, L1 = B1.T@B1 + B2@B2.T, but B2 is empty here)
L1 = B1.T @ B1
print_matrix("Edge Laplacian L1 = B1.T @ B1", L1)

# Let's compute L1 @ x1
L1_x1 = L1 @ x1

# As we derived, if x1 is harmonic, L1 @ x1 should be 0.
print("Option A check: Is x1 an eigenvector of L1?")
print("We compute the equation L1 @ x1 = lambda * x1")
print_matrix("L1 @ x1", L1_x1)

if np.allclose(L1_x1, 0):
    print("The result is the zero vector.")
    print("So, the equation holds for eigenvalue lambda = 0.")
    print("L1 @ x1 = 0 * x1")
    print("Final equation with numbers:")
    for i in range(num_edges):
        row_str = " + ".join([f"{L1[i, j]:.1f}*{x1[j]:.1f}" for j in range(num_edges)])
        print(f"(L1 @ x1)_{i}: {row_str} = {L1_x1[i]:.1f}")
    print("\nConclusion: x1 is an eigenvector of L1. Option A is true.")
else:
    print("x1 is not an eigenvector of L1. Option A is false.")
