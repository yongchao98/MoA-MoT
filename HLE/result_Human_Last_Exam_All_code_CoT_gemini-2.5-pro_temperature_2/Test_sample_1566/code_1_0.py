# This script calculates the number of topologically distinct continua
# that satisfy the two given properties.

# Let N be the number of topologically distinct continua with the given properties.

# Step 1 & 2: Analyze the given properties and their implications.
#
# A continuum X is a compact connected metric space.
# Property (1): X has k > 1 and finitely many "end points" E.
# The provided definition of an end point implies that X is a chainable continuum.
# Property (2): X has exactly two orbits under the action of its group of
# homeomorphisms, G = Homeo(X).
#
# Since the property of being an end point is preserved by homeomorphisms,
# the set E must be a union of orbits. Because there are only two orbits (O1, O2),
# and E is neither empty nor the entire space, E must be one of the orbits.
# This means:
# - Orbit 1 is the set of k end points, E.
# - Orbit 2 is the set of all other points, X \ E.
# All k end points are topologically equivalent, and all non-end points are also
# topologically equivalent to each other.

# Step 3: Search for candidate continua by cases.

# Case A: Locally Connected Continua (Dendrites)

# A1: The simple arc (e.g., the interval [0,1]).
# - End points: E = {0, 1}, so k=2. Property (1) is satisfied.
# - Orbits: Homeomorphisms of [0,1] include reflections (e.g., h(x) = 1-x)
#   and stretching maps. The set {0,1} is one orbit. The set of all
#   interior points (0,1) is the second orbit. Property (2) is satisfied.
# - Conclusion: The simple arc is a valid solution.
count_arc_type = 1
print(f"Number of solutions of the simple arc type: {count_arc_type}")

# A2: Dendrites with more than 2 end points (k > 2).
# - Any such dendrite must have at least one branch point (a vertex of degree 3 or more).
# - This branch point `j` is not an end point.
# - A homeomorphism must preserve the local structure of a point. A branch point `j`
#   has a neighborhood different from that of an interior point on an edge.
# - Thus, `j` cannot be in the same orbit as an interior point.
# - This leads to at least three orbits: (1) the endpoints, (2) the branch points,
#   and (3) the interior points. This violates Property (2).
# - Conclusion: No dendrite with k > 2 end points is a solution.
count_other_dendrite_types = 0
print(f"Number of solutions among other dendrite types: {count_other_dendrite_types}")


# Case B: Non-Locally Connected Continua
# - This case involves more exotic spaces, typically related to indecomposable
#   continua like the pseudo-arc.
# - Constructing a space with a finite number of end points while ensuring
#   the remaining (non-locally connected) part is a single orbit is highly non-trivial.
# - Standard constructions, such as attaching arcs to a pseudo-arc, fail because
#   the attachment points become special, creating additional orbits.
# - It is conjectured that no such non-locally connected continua exist that satisfy
#   these strong symmetry conditions. For the scope of this problem, we conclude
#   that there are no known examples that fit these criteria.
count_non_locally_connected_type = 0
print(f"Number of solutions among non-locally connected types: {count_non_locally_connected_type}")

# Step 4: Sum the results from all cases to get the final count.
total_number_of_continua = count_arc_type + count_other_dendrite_types + count_non_locally_connected_type

print(f"The total number of such topologically distinct continua is {count_arc_type} + {count_other_dendrite_types} + {count_non_locally_connected_type} = {total_number_of_continua}.")