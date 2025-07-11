def solve_group_theory_question():
    """
    This function solves the user's question by referencing a known theorem
    in group theory and presenting the result.
    """

    # According to a theorem by Diananda and Yap (1969), the finite groups
    # containing a maximal by inclusion product-free set of size 2 are known.
    # The following is the complete list of these non-isomorphic groups.
    groups = [
        "Z_4 (cyclic group of order 4)",
        "Z_5 (cyclic group of order 5)",
        "Z_6 (cyclic group of order 6)",
        "Z_7 (cyclic group of order 7)",
        "Z_9 (cyclic group of order 9)",
        "Z_2 x Z_2 (Klein four-group)",
        "Z_3 x Z_3 (direct product of two cyclic groups of order 3)",
        "S_3 (symmetric group of order 6)",
        "Q_8 (quaternion group of order 8)"
    ]

    print("The finite groups containing a maximal by inclusion product-free set of size 2 are:")
    for group in groups:
        print(f"- {group}")

    # Calculate the total count
    count = len(groups)
    
    # Create the equation string: 1 + 1 + ... + 1 = count
    sum_parts = ["1"] * count
    equation = " + ".join(sum_parts)
    
    print(f"\nIn total, there are {count} such groups.")
    print("The final count is derived from this list:")
    print(f"{equation} = {count}")

solve_group_theory_question()