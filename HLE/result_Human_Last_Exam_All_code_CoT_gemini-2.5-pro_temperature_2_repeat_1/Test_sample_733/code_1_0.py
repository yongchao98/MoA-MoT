def solve_group_theory_problem():
    """
    This function presents the solution to the question about finite groups
    with maximal by inclusion product-free sets of size 2.
    
    The solution is based on a known classification result from finite group theory.
    """
    
    groups = [
        "C₄ (Cyclic group of order 4)",
        "V₄ (Z₂ x Z₂, the Klein four-group)",
        "C₅ (Cyclic group of order 5)",
        "C₆ (Cyclic group of order 6)",
        "S₃ (Symmetric group on 3 elements)",
        "D₈ (Dihedral group of order 8)",
        "Q₈ (Quaternion group)"
    ]
    
    number_of_groups = len(groups)
    
    print("The finite groups that contain a maximal by inclusion product-free set of size 2 are:")
    for group in groups:
        print(f"- {group}")
        
    print("\nTotal number of such groups is:")
    print(number_of_groups)

solve_group_theory_problem()