# Each of the five component spaces is topologically distinct from the others.
# An auto-homeomorphism of the total space must map each component space to itself.
# Therefore, the equivalence classes of the total space are the union of the
# equivalence classes of each component space.

# We calculate the number of equivalence classes for each space.
# A space is homogeneous if for any two points x, y, there exists an
# auto-homeomorphism sending x to y. A homogeneous space has exactly one
# equivalence class.

# 1. The torus is a homogeneous space.
num_classes_torus = 1

# 2. The sphere is a homogeneous space.
num_classes_sphere = 1

# 3. The real line is a homogeneous space.
num_classes_real_line = 1

# 4. Any finite discrete space is homogeneous, since any permutation of its
# points is an auto-homeomorphism.
num_classes_three_point_discrete = 1

# 5. The five-point discrete space is also homogeneous.
num_classes_five_point_discrete = 1

# The total number of equivalence classes is the sum of the number of classes
# in each disjoint component.
total_equivalence_classes = (
    num_classes_torus +
    num_classes_sphere +
    num_classes_real_line +
    num_classes_three_point_discrete +
    num_classes_five_point_discrete
)

# Print the final calculation and result.
print(f"The number of equivalence classes for each space are:")
print(f"Torus: {num_classes_torus}")
print(f"Sphere: {num_classes_sphere}")
print(f"Real line: {num_classes_real_line}")
print(f"Three-point discrete space: {num_classes_three_point_discrete}")
print(f"Five-point discrete space: {num_classes_five_point_discrete}")
print("\nThe total number of equivalence classes is the sum of these numbers:")
print(f"{num_classes_torus} + {num_classes_sphere} + {num_classes_real_line} + {num_classes_three_point_discrete} + {num_classes_five_point_discrete} = {total_equivalence_classes}")