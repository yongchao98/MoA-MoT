# The problem is to find the number of certain components of a fractal set F.
# The set F is defined by the equation F = union_{d in D} (F+d)/4,
# where D = {(0,0), (0,1), (0,2),(0,3),(3,0), (3,1), (3,2),(3,3)}.
# This set F is the attractor of an Iterated Function System (IFS).

# Step 1: Analyze the connectivity of the IFS.
# The connectivity of the attractor F is determined by the connectivity of a graph G.
# The vertices of G are the elements of D.
# An edge exists between d_i and d_j if the set images f_i(F) and f_j(F) intersect.
# The condition for intersection is that d_j - d_i can be written as x - y for some x, y in F.

# Step 2: Characterize the set F and the difference x - y.
# Since F is a subset of the unit square [0,1]^2, the difference vector x - y
# must lie in [-1, 1] x [-1, 1].

# Step 3: Analyze the differences d_j - d_i.
# Let d_i = (dx_i, dy_i) and d_j = (dx_j, dy_j).
# dx_k can be 0 or 3. dy_k can be 0, 1, 2, or 3.
# The difference vector is (dx_j - dx_i, dy_j - dy_i).
# The x-component of the difference can be -3, 0, or 3.
# The y-component of the difference can be -3, -2, -1, 0, 1, 2, or 3.

# Step 4: Apply the constraint from Step 2 to the differences in Step 3.
# For an intersection to occur, the difference vector must be in [-1, 1] x [-1, 1].
# This means:
# - dx_j - dx_i must be 0. This separates the vertices into two groups based on the x-coordinate.
# - dy_j - dy_i must be in {-1, 0, 1}.

# So, an edge exists between two distinct vectors d_i and d_j if and only if:
# dx_i = dx_j AND |dy_i - dy_j| = 1.

# Step 5: Construct the graph and count its connected components.
# The vertices are D = {(0,0), (0,1), (0,2), (0,3), (3,0), (3,1), (3,2), (3,3)}.
# Group 1 (x-coordinate is 0): V1 = {(0,0), (0,1), (0,2), (0,3)}
# Edges: ((0,0),(0,1)), ((0,1),(0,2)), ((0,2),(0,3)).
# This group forms a single connected component (a path).

# Group 2 (x-coordinate is 3): V2 = {(3,0), (3,1), (3,2), (3,3)}
# Edges: ((3,0),(3,1)), ((3,1),(3,2)), ((3,2),(3,3)).
# This group also forms a single connected component (a path).

# There are no edges between V1 and V2 because their x-coordinates differ.
# So, the graph has two connected components.

# Step 6: Conclude the number of components of F.
# According to the theory of self-similar sets, the number of connected components
# of the attractor F is equal to the number of connected components of this graph G.
# The properties of being "nondegenerate" and "locally connected" are assumed to hold for these components.
# Therefore, the number of such components is 2.

num_components = 2
print("The smallest possible number of components of F that are nondegenerate and locally connected is:")
print(num_components)