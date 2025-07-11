# The problem is to find the minimum total number of edges in a topologically nontrivial
# link of three components on the 3D integer lattice.

# Let the three components be K1, K2, and K3.
# The link is nontrivial if the components cannot all be separated.
# The minimum length is achieved by the simplest form of nontriviality:
# two components (K1, K2) are linked, and the third (K3) is separate.

# The minimal configuration for two linked components (a Hopf link) consists of:
# 1. A 2x2 square knot. Its path length is 2+2+2+2.
K1_edges = 8
# 2. A 1x1 square knot that passes through the hole of the 2x2 square.
#    Its path length is 1+1+1+1.
K2_edges = 4

# The combined length of this minimal Hopf link is 8 + 4 = 12.

# The third component (K3) can be a minimal unlinked knot.
# The smallest possible knot on the lattice is a 1x1 square.
K3_edges = 4

# The total minimum number of edges is the sum of the edges of these three components.
total_edges = K1_edges + K2_edges + K3_edges

# The final answer is the calculation showing how this minimum is achieved.
# The problem asks us to output each number in the final equation.
print(f"The minimum configuration is a Hopf link plus an unlinked component.")
print(f"Hopf link component 1 (e.g., a 2x2 square): {K1_edges} edges")
print(f"Hopf link component 2 (e.g., a 1x1 square): {K2_edges} edges")
print(f"Unlinked component 3 (e.g., a 1x1 square): {K3_edges} edges")
print(f"Total minimum edges = {K1_edges} + {K2_edges} + {K3_edges} = {total_edges}")