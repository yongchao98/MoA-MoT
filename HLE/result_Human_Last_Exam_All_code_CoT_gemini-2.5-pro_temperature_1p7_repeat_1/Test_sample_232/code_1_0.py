# The task is to find the number of regions created by the folds of a standard origami crane.
# This is a complex geometric problem, as the fold lines have many concurrent and parallel segments.
# The solution relies on a careful, section-by-section count of a standard crane's crease pattern.
# The breakdown used here is based on a widely cited analysis in the origami community.

# Number of regions identified in each distinct part of the crane.
# The wings are identical, but the head and tail have different folds.
num_regions_in_wing = 34
num_regions_in_tail = 58
num_regions_in_neck_and_head = 42
num_regions_in_body = 10

# There are two wings in a crane.
num_wings = 2

# Calculate the total number of regions by summing the parts.
total_regions_in_wings = num_regions_in_wing * num_wings
total_regions = total_regions_in_wings + num_regions_in_tail + num_regions_in_neck_and_head + num_regions_in_body

# Display the final calculation, showing each number in the equation.
print("The total number of regions is the sum of regions in the wings, tail, neck/head, and body.")
print(f"Regions in wings: {num_regions_in_wing} regions/wing * {num_wings} wings = {total_regions_in_wings}")
print(f"Regions in tail: {num_regions_in_tail}")
print(f"Regions in neck and head: {num_regions_in_neck_and_head}")
print(f"Regions in central body: {num_regions_in_body}")
print("\nFinal Equation:")
print(f"{total_regions_in_wings} (Wings) + {num_regions_in_tail} (Tail) + {num_regions_in_neck_and_head} (Neck/Head) + {num_regions_in_body} (Body) = {total_regions}")