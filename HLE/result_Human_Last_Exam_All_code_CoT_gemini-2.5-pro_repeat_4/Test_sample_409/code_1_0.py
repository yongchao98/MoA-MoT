def find_minimal_cohomology_degree():
    """
    This function analyzes the roles of different cohomology degrees in semi-abelian
    categories to find the minimal degree for classifying non-trivial extensions
    and obstructions.
    """

    # A dictionary mapping cohomology degrees to their primary role.
    cohomology_roles = {
        0: "H^0 classifies invariants or fixed points.",
        1: "H^1 classifies derivations and split extensions.",
        2: "H^2 classifies non-trivial extensions and primary obstructions.",
        3: "H^3 classifies more complex structures (e.g., crossed modules) and higher obstructions."
    }

    print("A summary of the roles of low-degree cohomology groups H^n(B, A):")
    for degree, role in cohomology_roles.items():
        print(f"Degree {degree}: {role}")

    # The question asks for the minimal degree where *non-trivial* extensions
    # AND obstructions become significant.
    # While H^1 is related to extensions, H^2 is the classical degree for
    # classifying the more complex, non-trivial extensions (like in group theory)
    # and is the primary locus for obstruction theory.
    minimal_degree = 2

    print("\nConclusion:")
    print("The minimal degree at which cohomology becomes significant for both non-trivial extensions and obstructions is 2.")

    # Printing the result in an equation format as requested.
    print("\nFinal Answer Equation:")
    print(f"Minimal Significant Cohomology Degree = {minimal_degree}")


if __name__ == "__main__":
    find_minimal_cohomology_degree()