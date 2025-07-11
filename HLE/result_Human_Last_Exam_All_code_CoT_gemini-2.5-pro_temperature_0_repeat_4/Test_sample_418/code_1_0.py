# This script calculates the number of equivalence classes for the given topological space.

# The space X is a disjoint union of 5 topologically distinct spaces.
# An auto-homeomorphism of X must map each component to itself.
# Therefore, the total number of equivalence classes is the sum of the
# number of equivalence classes in each component.

# Each of the five component spaces is homogeneous, meaning any point can be
# mapped to any other point within that component by an auto-homeomorphism.
# A homogeneous space has exactly one equivalence class.

# Number of equivalence classes for the Torus
num_classes_torus = 1

# Number of equivalence classes for the Sphere
num_classes_sphere = 1

# Number of equivalence classes for the Real line
num_classes_real_line = 1

# Number of equivalence classes for a three-point discrete space
num_classes_d3 = 1

# Number of equivalence classes for a five-point discrete space
num_classes_d5 = 1

# Total number of equivalence classes is the sum
total_classes = num_classes_torus + num_classes_sphere + num_classes_real_line + num_classes_d3 + num_classes_d5

# Print the equation showing the contribution from each space
print(f"The total number of equivalence classes is the sum of the classes from each component:")
print(f"{num_classes_torus} (Torus) + {num_classes_sphere} (Sphere) + {num_classes_real_line} (Real line) + {num_classes_d3} (3-pt discrete) + {num_classes_d5} (5-pt discrete) = {total_classes}")