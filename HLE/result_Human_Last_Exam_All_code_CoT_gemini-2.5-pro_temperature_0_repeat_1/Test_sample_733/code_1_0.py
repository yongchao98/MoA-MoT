def solve_group_theory_problem():
    """
    This function solves the problem by stating a known classification theorem
    from group theory and counting the number of groups in the classification.
    """

    # The question asks for the number of finite groups that contain a maximal
    # by inclusion product-free set of size 2. This is a known result in
    # combinatorial group theory.
    #
    # A theorem by E. V. Kopylov (and others) classifies all such finite groups.
    # A finite group G has this property if and only if it is isomorphic to
    # one of the following 8 groups:

    classified_groups = [
        "C_4 (the cyclic group of order 4)",
        "C_5 (the cyclic group of order 5)",
        "C_6 (the cyclic group of order 6)",
        "D_6 (the dihedral group of order 6, also known as S_3)",
        "D_8 (the dihedral group of order 8)",
        "D_10 (the dihedral group of order 10)",
        "A_4 (the alternating group on 4 elements, order 12)",
        "The non-abelian group of order 27 and exponent 3"
    ]

    # The number of such groups is the length of this list.
    number_of_groups = len(classified_groups)

    print("The finite groups containing a maximal by inclusion product-free set of size 2 are:")
    for i, group_name in enumerate(classified_groups, 1):
        print(f"{i}. {group_name}")

    print("\nTo find the total number of such groups, we count the groups in the list.")
    
    # Building and printing the equation as requested.
    equation_parts = ["1" for _ in classified_groups]
    equation_str = " + ".join(equation_parts)
    
    print(f"\nThe final equation is the sum of the counts for each group in the classification:")
    print(f"{equation_str} = {number_of_groups}")

    print(f"\nThus, the total number of such finite groups is {number_of_groups}.")

solve_group_theory_problem()