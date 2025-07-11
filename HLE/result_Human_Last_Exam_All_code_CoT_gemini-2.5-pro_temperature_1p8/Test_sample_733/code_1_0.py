def solve_group_theory_question():
    """
    This function provides the solution to the question based on a known classification
    theorem from group theory.
    """
    # The classification theorem by R. M. Green states that a finite group G has a
    # maximal by inclusion product-free set of size 2 if and only if it is
    # isomorphic to one of the following five groups.
    # We list their standard names.
    
    groups = [
        "C_4 (the cyclic group of order 4)",
        "C_5 (the cyclic group of order 5)",
        "S_3 (the symmetric group on 3 elements, also known as D_6, the dihedral group of order 6)",
        "D_8 (the dihedral group of order 8, symmetries of a square)",
        "D_10 (the dihedral group of order 10, symmetries of a pentagon)"
    ]

    print("The finite groups containing maximal by inclusion product-free sets of size 2 are:")
    for group in groups:
        print(f"- {group}")
        
    count = len(groups)
    
    print("\nThe number of such non-isomorphic groups is:")
    print(count)

solve_group_theory_question()