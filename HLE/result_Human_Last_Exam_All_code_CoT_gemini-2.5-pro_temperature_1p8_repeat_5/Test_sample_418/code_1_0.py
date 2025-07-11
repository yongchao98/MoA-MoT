# Number of equivalence classes for each component space.
# Each of these spaces is homogeneous, meaning all points within it are equivalent.
# Therefore, each space contributes exactly one equivalence class.

# Number of equivalence classes for the Torus
torus_classes = 1

# Number of equivalence classes for the Sphere
sphere_classes = 1

# Number of equivalence classes for the Real line
real_line_classes = 1

# Number of equivalence classes for a three-point discrete space
three_point_discrete_classes = 1

# Number of equivalence classes for a five-point discrete space
five_point_discrete_classes = 1

# The total number of equivalence classes is the sum of the classes from each
# component, as any auto-homeomorphism of the total space must map each
# component to itself.
total_classes = (torus_classes +
                 sphere_classes +
                 real_line_classes +
                 three_point_discrete_classes +
                 five_point_discrete_classes)

# Print the final equation and the result
print("The total number of equivalence classes is the sum of the classes from each component.")
print(f"The calculation is: {torus_classes} (Torus) + {sphere_classes} (Sphere) + {real_line_classes} (Real Line) + {three_point_discrete_classes} (3-point discrete) + {five_point_discrete_classes} (5-point discrete) = {total_classes}")
