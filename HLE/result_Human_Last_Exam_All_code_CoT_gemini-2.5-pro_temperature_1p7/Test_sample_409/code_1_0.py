def explain_cohomology(degree):
    """
    Provides a brief explanation for the role of the n-th cohomology group
    H^n(B, A) in the context of extensions in semi-abelian categories.
    """
    role = ""
    if degree == 0:
        role = "H^0(B, A): Classifies the invariants of the B-module A. Not directly related to extensions."
    elif degree == 1:
        role = "H^1(B, A): Classifies derivations and splittings of semi-direct products. Does not classify general, non-split extensions."
    elif degree == 2:
        role = "H^2(B, A): Classifies all extensions of B by A. A non-zero element represents a non-trivial (non-split) extension, acting as an obstruction to splitting."
    elif degree == 3:
        role = "H^3(B, A): Classifies higher-order structures and obstructions, such as those related to crossed modules."
    else:
        role = f"H^{degree}(B, A) deals with even more abstract algebraic structures."

    print(role)

def solve_cohomology_question():
    """
    Analyzes the roles of low-degree cohomology groups to find the minimal
    degree for classifying non-trivial extensions and obstructions.
    """
    print("Analyzing the role of cohomology groups H^n(B, A) by degree:\n")
    for i in range(4):
        explain_cohomology(i)

    print("\n--- Conclusion ---")
    print("The analysis shows that degree 2 is the minimal degree where non-trivial extensions and obstructions to splitting them are fundamentally classified.")

    # Per the instructions, we print the final equation showing the numbers involved.
    print("\nThis fundamental classification result is captured by the following relationship:")
    eq_text_1 = "The set of equivalence classes of extensions of B by A is in bijection with"
    eq_part_1 = "H"
    eq_part_2 = 2
    eq_part_3 = "(B, A)"
    print(f"{eq_text_1} {eq_part_1}^{eq_part_2}{eq_part_3}")
    print(f"The significant number in this central equation of extension theory is: {eq_part_2}")

# Execute the explanatory function
solve_cohomology_question()