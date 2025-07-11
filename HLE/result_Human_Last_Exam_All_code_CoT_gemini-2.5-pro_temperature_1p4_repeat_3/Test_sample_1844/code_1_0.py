# 1. Define the number of edges for the outer perimeter.
# A square has 4 sides. Cutting the 4 corners each adds an edge,
# and truncates the original sides, but the total count of edges for the
# outer perimeter becomes 8.
outer_perimeter_edges = 8

# 2. Define the number of edges for the central hole.
# The cut at the centermost corner unfolds to a 4-sided hole.
central_hole_edges = 4

# 3. Define the number of edges for the holes on the horizontal axis.
# The cut at another corner creates two 4-sided holes along the horizontal axis.
h_axis_holes_edges = 2 * 4

# 4. Define the number of edges for the holes on the vertical axis.
# The cut at the final corner creates two 4-sided holes along the vertical axis.
v_axis_holes_edges = 2 * 4

# 5. Calculate the total number of edges.
total_edges = outer_perimeter_edges + central_hole_edges + h_axis_holes_edges + v_axis_holes_edges

# Print the breakdown of the calculation as requested.
print("The total number of edges is the sum of:")
print(f"- Edges of the outer perimeter: {outer_perimeter_edges}")
print(f"- Edges of the central hole: {central_hole_edges}")
print(f"- Edges of holes on the horizontal axis: {h_axis_holes_edges}")
print(f"- Edges of holes on the vertical axis: {v_axis_holes_edges}")
print("-" * 20)
# The final output prints each number in the equation.
print(f"Total Edges = {outer_perimeter_edges} + {central_hole_edges} + {h_axis_holes_edges} + {v_axis_holes_edges} = {total_edges}")
<<<28>>>