def explain_cohomology_degree():
    """
    Explains the roles of low-degree cohomology groups in semi-abelian categories
    to determine the minimal degree for non-trivial extensions and obstructions.
    """
    print("Step 1: Analyzing the meaning of low-degree cohomology groups H^n(B, A).")
    print("-" * 60)

    # --- Degree 0 ---
    degree_0 = 0
    print(f"For n = {degree_0}:")
    print("H^0(B, A) represents the 'invariants'. It consists of the elements of the module A that are 'fixed' by the action of B. This degree is not concerned with extensions or obstructions.")
    print("")

    # --- Degree 1 ---
    degree_1 = 1
    print(f"For n = {degree_1}:")
    print("H^1(B, A) classifies 'derivations' (or crossed homomorphisms). It is connected to the study of split extensions, but it does not classify the richer structure of general non-trivial extensions.")
    print("")

    # --- Degree 2 ---
    degree_2 = 2
    print(f"For n = {degree_2}:")
    print(f"H^2(B, A) is fundamentally important for both extensions and obstructions.")
    print("  - Non-trivial Extensions: This is the classical result from Eilenberg and MacLane. The elements of H^2(B, A) are in a one-to-one correspondence with the equivalence classes of (abelian) extensions of B by A. A non-zero element in H^2(B, A) corresponds to a non-split, i.e., non-trivial, extension.")
    print("  - Obstructions: The class in H^2(B, A) associated with an extension can be interpreted as the 'obstruction' to that extension being split. The extension is trivial (split) if and only if this obstruction vanishes (is the zero element of the group).")
    print("")

    # --- Degree 3 ---
    degree_3 = 3
    print(f"For n = {degree_3}:")
    print("H^3(B, A) typically appears as the obstruction to more complex problems, such as whether an outer action of B on another object can be lifted to an actual extension.")
    print("")

    # --- Conclusion ---
    print("Step 2: Conclusion.")
    print("-" * 60)
    print(f"The question asks for the minimal degree where both non-trivial extensions and obstructions are significant. While H^3 is famous for obstruction theory, H^2 is the first degree where both concepts appear:")
    print(f"1. It classifies non-trivial extensions.")
    print(f"2. It provides the obstruction to an extension being trivial (split).")
    print(f"\nTherefore, the minimal significant degree is {degree_2}.")

explain_cohomology_degree()