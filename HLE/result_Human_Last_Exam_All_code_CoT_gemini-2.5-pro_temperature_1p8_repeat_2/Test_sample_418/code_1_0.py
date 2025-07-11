# Number of equivalence classes for each component space.
# A space is homogeneous if its group of auto-homeomorphisms acts transitively on its points.
# If a space is homogeneous, it has exactly one equivalence class.

# The torus is a homogeneous space.
torus_classes = 1

# The sphere is a homogeneous space.
sphere_classes = 1

# The real line is a homogeneous space.
real_line_classes = 1

# A discrete space with n points is homogeneous because the symmetric group S_n of all permutations acts transitively on the points.
three_point_discrete_classes = 1
five_point_discrete_classes = 1

# The total number of equivalence classes in the disjoint union is the sum of the
# number of equivalence classes in each component, as homeomorphisms cannot map points
# between topologically distinct components.
total_classes = torus_classes + sphere_classes + real_line_classes + three_point_discrete_classes + five_point_discrete_classes

# Print the final equation showing how the total is calculated.
print(f"{torus_classes} + {sphere_classes} + {real_line_classes} + {three_point_discrete_classes} + {five_point_discrete_classes} = {total_classes}")
