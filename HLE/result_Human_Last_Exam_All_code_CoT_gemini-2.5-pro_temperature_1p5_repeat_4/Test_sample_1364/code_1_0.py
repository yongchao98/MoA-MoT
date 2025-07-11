# This script is designed to output the solution to the stated problem.
# The problem asks for the set of possible numbers of vertices a convex polyhedron P
# can have, given that it has quadrilateral projections on three planes in general position.

# Based on geometric analysis, we identify candidates by examining families of polyhedra:
# - V=4: A tetrahedron.
# - V=5: A pyramid with a quadrilateral base.
# - V=6: A bipyramid with a quadrilateral base (e.g., an octahedron).
# - V=8: A prism or frustum with a quadrilateral base (e.g., a cube).

# These specific constructions allow for three linearly independent projection
# directions that result in quadrilateral silhouettes. Other numbers of vertices,
# such as V=7 or V>8, generally fail this condition for simple polyhedral families.
# For example, a pentagonal prism (V=10) has a pentagonal projection from one direction.

# The derived set of possible numbers of vertices.
possible_n_vertices = {4, 5, 6, 8}

# Output the result. The instruction "output each number in the final equation"
# is interpreted as printing the elements of the derived solution set.
print("The set of possible numbers of vertices is:")
# We convert the set to a sorted list for ordered printing.
result = sorted(list(possible_n_vertices))
for number in result:
  print(number)