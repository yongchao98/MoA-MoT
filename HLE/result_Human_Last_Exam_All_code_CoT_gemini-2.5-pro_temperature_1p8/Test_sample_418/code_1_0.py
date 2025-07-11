# The problem is to find the number of equivalence classes on a space X,
# which is the disjoint union of 5 different topological spaces.
# The equivalence relation x ~ y means there is an auto-homeomorphism of X sending x to y.

# Step 1: Identify the components.
# The space X is a disjoint union of 5 components. Auto-homeomorphisms map components
# to homeomorphic components. As all 5 components are topologically distinct, any
# auto-homeomorphism must map each component to itself.
# This means we can count the equivalence classes on each component separately and add them up.

# Step 2: Count equivalence classes on each component.
# A space is 'homogeneous' if any point can be mapped to any other point by an auto-homeomorphism.
# If a space is homogeneous, it has exactly one equivalence class.

# The list of component spaces.
component_spaces = [
    "The torus",
    "The sphere",
    "The real line",
    "A three point discrete space",
    "A five point discrete space"
]

# All the given component spaces are homogeneous.
# - Torus: Homogeneous via translations.
# - Sphere: Homogeneous via rotations.
# - Real line: Homogeneous via translations.
# - n-point discrete space: Homogeneous via permutations (which are homeomorphisms).
# So, each component contributes exactly 1 equivalence class.
classes_per_component = [1, 1, 1, 1, 1]

# Step 3: Sum the number of classes.
total_number_of_classes = sum(classes_per_component)

# Step 4: Display the result as a detailed equation.
# The problem asks to output each number in the final equation.
equation_string = " + ".join(map(str, classes_per_component))
print(f"Number of classes from each component (Torus, Sphere, Real line, 3-point discrete, 5-point discrete):")
print(f"{equation_string}")
print("Total number of equivalence classes:")
print(f"{equation_string} = {total_number_of_classes}")