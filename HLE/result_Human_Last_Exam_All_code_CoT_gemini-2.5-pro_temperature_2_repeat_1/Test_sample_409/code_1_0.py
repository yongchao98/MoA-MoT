def explain_cohomology_degree():
    """
    Explains the significance of cohomology degrees and determines the minimal
    degree for non-trivial extensions and obstructions in semi-abelian categories.
    """

    # The roles of low-degree cohomology groups H^n(B, M)
    roles = {
        0: "Invariants of the action of B on M.",
        1: "Classification of derivations and split (trivial) extensions.",
        2: "Classification of non-trivial extensions. The non-vanishing itself is the first significant obstruction (to triviality).",
        3: "Classification of higher structures and obstructions to lifting structures classified by H^2."
    }

    # The question seeks the minimal degree for "non-trivial extensions" and "obstructions".
    # H^1 deals with extensions, but they are the 'trivial' kind.
    # H^2 is the first to deal with non-trivial extensions. These extensions represent an
    # obstruction to solving a lifting problem, so H^2 is the first degree where
    # both concepts are fundamentally significant.
    minimal_degree = 2

    print("Analyzing the roles of cohomology degrees H^n(B, M):")
    for degree, role in roles.items():
        print(f"  - H^{degree}: {role}")

    print("\nConclusion:")
    print("Degree 2 is the minimal degree where non-trivial extensions are classified,")
    print("and the existence of such extensions represents the first significant type of obstruction.")

    # Fulfilling the request to output the number from the final 'equation'.
    # Here, the equation is simply: Final Answer = 2.
    print("\nFinal Answer Equation:")
    print(f"Minimal Significant Degree = {minimal_degree}")
    # Outputting the number in the equation.
    print(f"The number in the final equation is: {minimal_degree}")

explain_cohomology_degree()