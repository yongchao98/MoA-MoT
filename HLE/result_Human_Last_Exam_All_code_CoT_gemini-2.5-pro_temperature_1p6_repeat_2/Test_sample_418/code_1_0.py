# Each of the 5 component spaces is homogeneous, meaning the group of
# auto-homeomorphisms acts transitively on it. This implies that for each
# component space, all of its points belong to a single equivalence class.
# Therefore, each component space contributes exactly one equivalence class to the total.

# Number of classes for the Torus
torus_classes = 1
# Number of classes for the Sphere
sphere_classes = 1
# Number of classes for the Real line
real_line_classes = 1
# Number of classes for the three-point discrete space
three_point_classes = 1
# Number of classes for the five-point discrete space
five_point_classes = 1

# The list of class counts for each space
class_counts = [
    torus_classes,
    sphere_classes,
    real_line_classes,
    three_point_classes,
    five_point_classes
]

# The total number of equivalence classes is the sum of the classes from each space.
total_classes = sum(class_counts)

# Create the equation string to display the calculation.
# We explicitly show each number in the sum.
equation_parts = [str(count) for count in class_counts]
equation = " + ".join(equation_parts) + f" = {total_classes}"

print("The number of equivalence classes for each space are as follows:")
print(f"Torus: {torus_classes}")
print(f"Sphere: {sphere_classes}")
print(f"Real line: {real_line_classes}")
print(f"Three-point discrete space: {three_point_classes}")
print(f"Five-point discrete space: {five_point_classes}")
print("\nThe total number of equivalence classes is the sum:")
print(equation)