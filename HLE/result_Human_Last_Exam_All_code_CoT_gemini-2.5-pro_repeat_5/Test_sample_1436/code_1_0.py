import math

# The symmetry factor (S) of a Feynman diagram is a combinatorial number
# that counts the number of ways the diagram's components can be permuted
# without changing the diagram's structure.

# For the specified diagram with 2 vertices and 4 propagators connecting them:
# 1. Symmetry from permuting vertices: There are 2 indistinguishable vertices.
#    Swapping them leaves the diagram unchanged. This gives a factor of 2!.
# 2. Symmetry from permuting propagators: There are 4 identical propagators
#    (lines) connecting the two vertices. Permuting these lines among
#    themselves leaves the diagram unchanged. This gives a factor of 4!.
#
# The total symmetry factor is the product of these individual factors.

# Number of indistinguishable vertices that can be permuted
num_permutable_vertices = 2

# Number of identical propagators connecting the same two vertices
num_permutable_propagators = 4

# Calculate the factorial for vertex permutations
vertex_symmetry_factor = math.factorial(num_permutable_vertices)

# Calculate the factorial for propagator permutations
propagator_symmetry_factor = math.factorial(num_permutable_propagators)

# The total symmetry factor is the product
total_symmetry_factor = vertex_symmetry_factor * propagator_symmetry_factor

# Print the step-by-step calculation
print("Calculating the symmetry factor for a diagram with 2 vertices connected by 4 propagators.")
print("The formula for the symmetry factor (S) is:")
print("S = (permutations of identical vertices) * (permutations of identical propagators)")
print(f"S = {num_permutable_vertices}! * {num_permutable_propagators}!")
print(f"S = {vertex_symmetry_factor} * {propagator_symmetry_factor}")
print(f"S = {total_symmetry_factor}")