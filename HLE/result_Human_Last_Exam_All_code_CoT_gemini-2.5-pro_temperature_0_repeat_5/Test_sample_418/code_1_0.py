# Step 1: Determine the number of equivalence classes for each component space.
# A space is homogeneous if for any two points x, y, there is an auto-homeomorphism
# mapping x to y. Homogeneous spaces have exactly one equivalence class.

# The torus is homogeneous.
num_classes_torus = 1

# The sphere is homogeneous.
num_classes_sphere = 1

# The real line is homogeneous.
num_classes_real_line = 1

# Any discrete space is homogeneous because any permutation of its points is a homeomorphism.
num_classes_discrete_3 = 1
num_classes_discrete_5 = 1

# Step 2: The total number of equivalence classes is the sum of the classes
# from each topologically distinct component.
total_classes = (num_classes_torus +
                 num_classes_sphere +
                 num_classes_real_line +
                 num_classes_discrete_3 +
                 num_classes_discrete_5)

# Step 3: Print the final equation showing each number.
print(f"{num_classes_torus} + {num_classes_sphere} + {num_classes_real_line} + {num_classes_discrete_3} + {num_classes_discrete_5} = {total_classes}")