# Step 1: Define the number of equivalence classes for each component space.
# As explained in the reasoning, each of the five component spaces is homogeneous,
# meaning that all its points belong to a single equivalence class.

# The Torus is a homogeneous space.
num_classes_torus = 1

# The Sphere is a homogeneous space.
num_classes_sphere = 1

# The Real Line is a homogeneous space.
num_classes_real_line = 1

# A three-point discrete space is homogeneous.
num_classes_discrete_3 = 1

# A five-point discrete space is homogeneous.
num_classes_discrete_5 = 1

# Step 2: Sum the number of classes from each component.
# An auto-homeomorphism of the total space must map components to themselves,
# as they are all topologically distinct. Thus, the equivalence classes
# of the total space are the union of the equivalence classes of the components.
all_classes = [
    num_classes_torus,
    num_classes_sphere,
    num_classes_real_line,
    num_classes_discrete_3,
    num_classes_discrete_5
]

total_classes = sum(all_classes)

# Step 3: Print the final calculation, showing each term in the sum.
equation_str_parts = [str(n) for n in all_classes]
equation_str = " + ".join(equation_str_parts)

print(f"The total number of equivalence classes is the sum of the classes from each of the 5 component spaces.")
print(f"Total = {equation_str} = {total_classes}")
