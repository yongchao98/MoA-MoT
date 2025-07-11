def explain_cohomology_degrees():
    """
    Explains the significance of low-degree cohomology groups in semi-abelian categories
    to determine the minimal degree for non-trivial extensions and obstructions.
    """
    print("Analysis of Cohomology Degrees H^n(B, A)")
    print("=" * 40)
    print("Context: B is an object in a semi-abelian category, A is a B-module.\n")

    # Degree 0
    n = 0
    print(f"--- Degree n = {n} ---")
    print(f"H^{n}(B, A): Classifies B-invariant elements of A.")
    print("Significance: Represents stable elements or fixed points. It does not describe extensions or obstructions.\n")

    # Degree 1
    n = 1
    print(f"--- Degree n = {n} ---")
    print(f"H^{n}(B, A): Classifies split extensions of B by A (or equivalently, derivations).")
    print("Significance: This is the first degree where the concept of 'extensions' appears. However, it only captures the simplest, 'trivial' class of extensions that can be split.\n")

    # Degree 2
    n = 2
    print(f"--- Degree n = {n} ---")
    print(f"H^{n}(B, A): Classifies the general (potentially non-split) extensions of B by A.")
    print("This is represented by the set of short exact sequences:")
    print("  0 -> A -> E -> B -> 0")
    print("Significance: This is the minimal degree where truly non-trivial extensions are classified.")
    print("A non-zero element in H^2(B, A) corresponds to a non-split extension.")
    print("This non-triviality is precisely the OBSTRUCTION to the extension being a simpler split extension (semi-direct product).")
    print("Thus, at degree 2, both non-trivial extensions and obstructions become fundamentally significant.\n")

    # Degree 3
    n = 3
    print(f"--- Degree n = {n} ---")
    print(f"H^{n}(B, A): Classifies more complex structures and is the home of classical 'obstruction theory'.")
    print("Significance: It often measures the obstruction to extending other algebraic structures (like crossed modules). While crucial, the most fundamental obstruction (to splitting an extension) already appears in degree 2.\n")

    print("=" * 40)
    print("Conclusion:")
    print("Degree 2 is the minimal degree where the classification of non-trivial extensions begins,")
    print("and the cohomology classes themselves represent the obstructions to these extensions being split.")

    final_answer_degree = 2
    print("\nFinal Equation:")
    print(f"Minimal Significant Degree = {final_answer_degree}")


if __name__ == '__main__':
    explain_cohomology_degrees()