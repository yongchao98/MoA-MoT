# Each of the five spaces in the disjoint union is topologically distinct.
# An auto-homeomorphism of the total space must map each of these
# component spaces to itself. Therefore, the total number of equivalence
# classes is the sum of the number of classes in each component space.

# For a homogeneous space, all points belong to a single equivalence class.
# The number of equivalence classes is 1.

# The torus is a homogeneous space.
torus_classes = 1

# The sphere is a homogeneous space.
sphere_classes = 1

# The real line is a homogeneous space.
real_line_classes = 1

# A discrete space is homogeneous because any permutation of its points is
# a homeomorphism, and the permutation group acts transitively.
three_point_discrete_classes = 1
five_point_discrete_classes = 1

# The total number of equivalence classes is the sum.
total_classes = (torus_classes +
                 sphere_classes +
                 real_line_classes +
                 three_point_discrete_classes +
                 five_point_discrete_classes)

# Print the final equation showing the contribution from each space.
print(f"{torus_classes} (Torus) + {sphere_classes} (Sphere) + {real_line_classes} (Real Line) + {three_point_discrete_classes} (3-Point Discrete) + {five_point_discrete_classes} (5-Point Discrete) = {total_classes}")