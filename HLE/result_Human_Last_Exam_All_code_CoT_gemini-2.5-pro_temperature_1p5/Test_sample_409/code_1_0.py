def find_minimal_cohomology_degree():
    """
    This function analyzes the role of low-degree cohomology groups in semi-abelian
    categories to determine the minimal degree where non-trivial extensions and
    obstructions become significant.
    """
    
    print("Step-by-step analysis of cohomology degrees H^n(B, A):")
    
    # The prompt requires outputting numbers from the "final equation".
    # We will interpret this by explaining the role of each relevant degree number.

    # Degree 0
    degree_0 = 0
    role_0 = "Describes the B-invariant elements in A. It is not primarily concerned with extensions or obstructions."
    print(f"\nFor n = {degree_0}:")
    print(f"H^{degree_0} -> {role_0}")

    # Degree 1
    degree_1 = 1
    role_1 = "Classifies split extensions of B by A. While it concerns extensions, it does not capture the general case or significant obstruction theory."
    print(f"\nFor n = {degree_1}:")
    print(f"H^{degree_1} -> {role_1}")

    # Degree 2
    degree_2 = 2
    role_2 = "Classifies general (non-split) extensions of B by A. Its elements can also be seen as the obstruction to an extension being split. This is the first degree where both concepts are fully significant."
    print(f"\nFor n = {degree_2}:")
    print(f"H^{degree_2} -> {role_2}")

    # Degree 3
    degree_3 = 3
    role_3 = "Classifies higher-order structures and obstructions, such as the obstruction to the existence of an extension with specific properties."
    print(f"\nFor n = {degree_3}:")
    print(f"H^{degree_3} -> {role_3}")
    
    final_conclusion = 2
    print("\n-----------------------------------------------------------------")
    print(f"Conclusion: The minimal degree at which both non-trivial extensions")
    print(f"and obstructions are significant is {final_conclusion}.")
    print("-----------------------------------------------------------------")

find_minimal_cohomology_degree()