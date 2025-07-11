# Define a dictionary mapping cohomology degrees to their primary significance.
# This summarizes the standard interpretation in group cohomology, which is generalized
# by cohomology in semi-abelian categories.
cohomology_significance = {
    0: "Invariants (fixed points/subobjects)",
    1: "Derivations and Split Extensions",
    2: "Classification of general non-split extensions and obstructions to splitting",
    3: "Obstructions to extensions of crossed modules",
    4: "Higher obstructions and classifications"
}

# The problem asks for the minimal degree where both "non-trivial extensions"
# and "obstructions" become significant.
target_concepts = ["extensions", "obstructions"]
minimal_degree = -1

# We iterate through the sorted keys (the degrees) to find the minimum one
# whose description contains both key concepts.
for degree in sorted(cohomology_significance.keys()):
    description = cohomology_significance[degree].lower() # use lower for case-insensitive match
    if all(concept in description for concept in target_concepts):
        minimal_degree = degree
        break

print(f"Finding the minimal degree for 'extensions' and 'obstructions'...")
print(f"Degree {minimal_degree} has the description: '{cohomology_significance[minimal_degree]}'")
print(f"This is the minimal degree where both concepts are central.")
print(f"\nFinal Answer: {minimal_degree}")

# The question also mentions an "equation". The canonical structure classified by H^2
# is the short exact sequence representing an extension.
if minimal_degree == 2:
    print("\nThe key structure is the short exact sequence of an extension:")
    equation = "0 -> M -> E -> B -> 1"
    print(equation)
    print("\nThe numbers from this symbolic equation are:")
    # Printing each number from the equation string "0 -> ... -> 1"
    print(0)
    print(1)