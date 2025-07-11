# The hypothetical axiom states that through a point not on a line,
# there are exactly this many parallel lines.
num_parallels_per_vertex = 3

# We have a triangle ABC. We draw parallel lines from each vertex
# to the opposite side.

# Number of parallel lines drawn through vertex A (parallel to side BC)
lines_from_A = num_parallels_per_vertex

# Number of parallel lines drawn through vertex B (parallel to side AC)
lines_from_B = num_parallels_per_vertex

# Number of parallel lines drawn through vertex C (parallel to side AB)
lines_from_C = num_parallels_per_vertex

# Calculate the intersections between these sets of new lines.
# We assume that lines parallel to intersecting lines will also intersect.

# Intersections between lines from Group A and lines from Group B
intersections_A_vs_B = lines_from_A * lines_from_B

# Intersections between lines from Group A and lines from Group C
intersections_A_vs_C = lines_from_A * lines_from_C

# Intersections between lines from Group B and lines from Group C
intersections_B_vs_C = lines_from_B * lines_from_C

# The total number of new intersection points is the sum of these.
# Intersections within a group (e.g., two lines from Group A) or with the
# original triangle's sides only occur at the original vertices A, B, or C,
# which are excluded by the problem statement.
total_intersections = intersections_A_vs_B + intersections_A_vs_C + intersections_B_vs_C

print("Step-by-step calculation of new intersection points:")
print(f"Number of parallel lines per vertex: {num_parallels_per_vertex}")
print(f"Intersections between lines from A and B: {lines_from_A} * {lines_from_B} = {intersections_A_vs_B}")
print(f"Intersections between lines from A and C: {lines_from_A} * {lines_from_C} = {intersections_A_vs_C}")
print(f"Intersections between lines from B and C: {lines_from_B} * {lines_from_C} = {intersections_B_vs_C}")
print(f"Total new intersection points = {intersections_A_vs_B} + {intersections_A_vs_C} + {intersections_B_vs_C} = {total_intersections}")

print("\nFinal Answer:")
print(f"The total number of distinct points of intersection is {total_intersections}.")
<<<27>>>