# The set of possible numbers of vertices for a convex polyhedron P
# such that it has quadrilateral projections on three planes in general position.

# Based on geometric constructions, the possible values are:
# n=4: A tetrahedron.
# n=5: A square pyramid.
# n=6: An octahedron.
# n=7: An octahedron with a shallow cap on one face.
# n=8: A cube or a frustum of a square pyramid.

# Numbers of vertices greater than 8 are conjectured to be impossible.

possible_n_values = [4, 5, 6, 7, 8]

print("The set of possible numbers of vertices is:")
# The prompt requires outputting each number, so we iterate and print.
for n in possible_n_values:
    print(n)
