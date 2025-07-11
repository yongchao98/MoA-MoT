def find_minimal_cohomology_degree():
    """
    Analyzes the significance of low-degree cohomology in semi-abelian categories.

    In semi-abelian cohomology theory for B-modules:
    - H^0(B, A): Corresponds to invariants.
    - H^1(B, A): Classifies split extensions (a specific type of extension).
    - H^2(B, A): Classifies general extensions and introduces the first significant
                 notion of obstruction (e.g., the obstruction to an extension being split).
    - H^3(B, A): Deals with higher-order obstructions.

    The question asks for the minimal degree where both non-trivial extensions
    and obstructions become significant. While H^1 deals with extensions, H^2 is
    the first degree to handle general extensions and the associated obstructions
    in a fundamental way.
    """
    # The minimal degree where general extensions and obstructions are significant.
    minimal_degree = 2

    print("In semi-abelian categories, the roles of the first few cohomology degrees are:")
    print("H^0: Invariants")
    print("H^1: Derivations and Split Extensions")
    print("H^2: General Extensions and Obstructions")
    print("H^3: Higher Obstructions")
    print("\nThe minimal degree where both non-trivial extensions (in the general sense) and obstructions become significant is 2.")

    # The prompt requires printing the number from the final "equation".
    print("\nFinal Equation:")
    print(f"minimal_significant_degree = {minimal_degree}")

find_minimal_cohomology_degree()