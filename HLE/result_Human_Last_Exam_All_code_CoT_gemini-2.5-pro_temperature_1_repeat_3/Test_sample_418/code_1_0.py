# The problem asks for the number of equivalence classes for a topological space X,
# where x ~ y if there's an auto-homeomorphism of X mapping x to y.
# X is the disjoint union of: Torus, Sphere, Real Line, 3-point discrete space,
# and 5-point discrete space.

# Any auto-homeomorphism of X must map its connected components (the five spaces)
# to themselves, as they are not homeomorphic to each other.
# Therefore, the equivalence classes of X are the equivalence classes within each component.
# We need to find the number of equivalence classes for each component and sum them up.

# A space is "homogeneous" if for any two points x, y, there is an auto-homeomorphism
# mapping x to y. A homogeneous space has exactly one equivalence class.

# Number of equivalence classes in the Torus:
# The torus is homogeneous (e.g., via translations).
num_classes_torus = 1

# Number of equivalence classes in the Sphere:
# The sphere is homogeneous (e.g., via rotations).
num_classes_sphere = 1

# Number of equivalence classes in the Real Line:
# The real line is homogeneous (e.g., via translations).
num_classes_real_line = 1

# Number of equivalence classes in a three-point discrete space:
# Any permutation of points is a homeomorphism in a discrete space.
# Thus, it's homogeneous.
num_classes_discrete_3 = 1

# Number of equivalence classes in a five-point discrete space:
# Similarly, any permutation of points is a homeomorphism.
# Thus, it's also homogeneous.
num_classes_discrete_5 = 1

# The total number of equivalence classes is the sum of the classes from each component.
total_classes = (num_classes_torus +
                 num_classes_sphere +
                 num_classes_real_line +
                 num_classes_discrete_3 +
                 num_classes_discrete_5)

# Print the final equation as requested.
print(f"{num_classes_torus} + {num_classes_sphere} + {num_classes_real_line} + {num_classes_discrete_3} + {num_classes_discrete_5} = {total_classes}")
