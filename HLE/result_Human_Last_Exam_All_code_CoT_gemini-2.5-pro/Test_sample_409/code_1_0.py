def solve_cohomology_question():
    """
    Determines the minimal cohomology degree for significant non-trivial extensions
    and obstructions in semi-abelian categories.
    """

    # In homological algebra and its generalizations to semi-abelian categories,
    # low-degree cohomology groups have specific interpretations.
    cohomology_interpretations = {
        0: "H^0: Invariants. Does not describe extensions or obstructions.",
        1: "H^1: Derivations and split extensions. Does not describe general, non-split extensions.",
        2: "H^2: Classification of all extensions. A non-zero class is an obstruction to an extension being split.",
        3: "H^3: Obstruction theory (e.g., to the existence of an extension with certain properties)."
    }

    # The user asks for the MINIMAL degree where "non-trivial extensions" AND "obstructions"
    # become significant.

    # H^2 is the first degree that classifies non-trivial (i.e., non-split) extensions.
    # The existence of a non-zero class in H^2 for a given B and M is an
    # obstruction for every extension of B by M to be a split one (a semi-direct product).
    # Therefore, degree 2 is the minimal degree where both concepts are fundamentally significant.

    minimal_significant_degree = 2

    print(f"The roles of the first few cohomology degrees are:")
    for degree, description in cohomology_interpretations.items():
        print(f"  - {description}")

    print(f"\nThe minimal degree for classifying non-trivial extensions is {minimal_significant_degree}.")
    print(f"This degree, {minimal_significant_degree}, is also the first where the concept of obstruction (to an extension being split) becomes central.")
    print(f"Therefore, the final answer is {minimal_significant_degree}.")

solve_cohomology_question()