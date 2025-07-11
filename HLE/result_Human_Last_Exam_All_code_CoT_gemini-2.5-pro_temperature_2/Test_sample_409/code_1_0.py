def solve_cohomology_question():
    """
    This function determines the minimal significant cohomology degree
    for non-trivial extensions and obstructions in semi-abelian categories.
    """
    
    # In semi-abelian cohomology (generalizing group cohomology):
    # H^1 classifies split extensions.
    # H^2 classifies non-split extensions.
    # Therefore, the first degree where non-trivial (non-split) extensions are classified is 2.
    degree_for_nontrivial_extensions = 2
    
    # An element in H^2 that corresponds to a non-split extension can be viewed
    # as the obstruction to that extension being split.
    # More complex obstructions appear in H^3, but the most fundamental one starts here.
    # Thus, the minimal degree where obstructions become significant is also 2.
    degree_for_obstructions = 2
    
    # The minimal degree for either phenomenon is the minimum of these.
    minimal_significant_degree = min(degree_for_nontrivial_extensions, degree_for_obstructions)
    
    print("In cohomology theory, the significance of low-degree groups is as follows:")
    print("H^0: Invariants")
    print("H^1: Derivations / Split Extensions")
    print("H^2: Classification of Non-Split Extensions / Obstruction to splitting")
    print("H^3: Higher Obstruction Theory\n")
    
    print("The minimal degree where non-trivial extensions are classified is 2.")
    print("The minimal degree where significant obstructions appear is 2.")
    
    # Formulating this as an equation:
    # Minimal Degree = min(Degree(Extensions), Degree(Obstructions))
    equation_str = f"{minimal_significant_degree} = min({degree_for_nontrivial_extensions}, {degree_for_obstructions})"
    
    print("\nSymbolic Equation:")
    print(equation_str)
    
    print("\nAs requested, here is each number from the final equation:")
    print(f"Final Result: {minimal_significant_degree}")
    print(f"Number for Extensions: {degree_for_nontrivial_extensions}")
    print(f"Number for Obstructions: {degree_for_obstructions}")

solve_cohomology_question()