# Number of equivalence classes for the torus.
# The torus is a homogeneous space, so all points belong to one class.
torus_classes = 1

# Number of equivalence classes for the sphere.
# The sphere is a homogeneous space.
sphere_classes = 1

# Number of equivalence classes for the real line.
# The real line is a homogeneous space.
real_line_classes = 1

# Number of equivalence classes for a three point discrete space.
# Any n-point discrete space is homogeneous as any permutation of points is a homeomorphism.
three_point_classes = 1

# Number of equivalence classes for a five point discrete space.
# This space is also homogeneous.
five_point_classes = 1

# The total space is a disjoint union of non-homeomorphic components.
# Therefore, the total number of equivalence classes is the sum of the classes from each component.
total_classes = torus_classes + sphere_classes + real_line_classes + three_point_classes + five_point_classes

# Print the final equation and the result as requested.
print(f"{torus_classes} + {sphere_classes} + {real_line_classes} + {three_point_classes} + {five_point_classes} = {total_classes}")
