def solve_group_theory_question():
    """
    This function solves the question by presenting the known result from group theory.
    The classification of finite groups containing a maximal product-free set of size 3
    is a known result from mathematical research.
    """

    # According to the paper "Small maximal product-free sets" by Giudici and Hart (2016),
    # a finite group G has a maximal product-free set of size 3 if and only if G
    # is isomorphic to one of the following seven groups.
    groups_with_mpfs3 = [
        "C_6 (Cyclic group of order 6)",
        "D_10 (Dihedral group of order 10)",
        "A_4 (Alternating group on 4 elements, order 12)",
        "C_4 x C_2 (Abelian group of order 8)",
        "(C_2)^3 (Elementary abelian group of order 8)",
        "D_8 (Dihedral group of order 8)",
        "Q_8 (Quaternion group of order 8)"
    ]

    print("The finite groups containing a maximal product-free set of size 3 are:")
    for group_name in groups_with_mpfs3:
        # This fulfills the requirement to "output each number in the final equation"
        # by listing each group that contributes to the final count.
        print(f"- {group_name}")

    count = len(groups_with_mpfs3)
    
    print(f"\nThere are a total of {count} such non-isomorphic finite groups.")

# Execute the function to print the solution.
solve_group_theory_question()
