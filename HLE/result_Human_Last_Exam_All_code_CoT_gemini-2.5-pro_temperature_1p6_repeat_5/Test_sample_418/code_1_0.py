# The space X is a disjoint union of five component spaces.
# The total number of equivalence classes is the sum of the number of classes
# in each component, since no component is homeomorphic to another.

# Number of equivalence classes for the Torus.
# The torus is a homogeneous space, so all points are equivalent.
classes_torus = 1

# Number of equivalence classes for the Sphere.
# The sphere is a homogeneous space, so all points are equivalent.
classes_sphere = 1

# Number of equivalence classes for the Real Line.
# The real line is a homogeneous space, so all points are equivalent.
classes_real_line = 1

# Number of equivalence classes for a three-point discrete space.
# Any permutation of points in a discrete space is a homeomorphism,
# so the space is homogeneous.
classes_3pt_discrete = 1

# Number of equivalence classes for a five-point discrete space.
# This space is also homogeneous for the same reason.
classes_5pt_discrete = 1

# Calculate the total number of equivalence classes by summing them up.
total_classes = classes_torus + classes_sphere + classes_real_line + classes_3pt_discrete + classes_5pt_discrete

# Print the final equation with each number.
print(f"The total number of equivalence classes is the sum of the classes from each of the 5 spaces.")
print(f"Total Classes = {classes_torus} + {classes_sphere} + {classes_real_line} + {classes_3pt_discrete} + {classes_5pt_discrete} = {total_classes}")

<<<5>>>