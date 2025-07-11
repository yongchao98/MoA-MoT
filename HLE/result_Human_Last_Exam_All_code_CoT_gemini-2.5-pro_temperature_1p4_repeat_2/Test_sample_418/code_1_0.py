# Number of equivalence classes for each component space
# Each of the spaces is homogeneous, meaning for any two points x, y in the space,
# there is an auto-homeomorphism that sends x to y.
# This implies that each space consists of a single equivalence class.

torus_classes = 1
sphere_classes = 1
real_line_classes = 1
three_point_discrete_classes = 1
five_point_discrete_classes = 1

# The total number of equivalence classes is the sum of the number of classes
# from each of the disjoint component spaces.
total_classes = (torus_classes + 
                 sphere_classes + 
                 real_line_classes + 
                 three_point_discrete_classes + 
                 five_point_discrete_classes)

# Print the final equation
print(f"{torus_classes} + {sphere_classes} + {real_line_classes} + {three_point_discrete_classes} + {five_point_discrete_classes} = {total_classes}")