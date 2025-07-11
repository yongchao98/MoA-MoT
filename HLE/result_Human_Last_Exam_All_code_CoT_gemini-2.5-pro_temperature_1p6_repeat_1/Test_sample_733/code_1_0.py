def solve_group_problem():
    """
    This function presents the solution to the question about finite groups
    with maximal by inclusion product-free sets of size 2.
    """
    # The classification of finite groups containing a maximal product-free set
    # of size 2 is a known result in mathematics.
    groups = [
        "The cyclic group of order 4 (Z_4)",
        "The Klein four-group (Z_2 x Z_2)",
        "The cyclic group of order 5 (Z_5)",
        "The cyclic group of order 6 (Z_6)",
        "The symmetric group on 3 elements (S_3 or D_3)"
    ]

    print("The finite groups that contain a maximal by inclusion product-free set of size 2 are:")
    for group in groups:
        print(f"- {group}")

    count = len(groups)
    print(f"\nIn total, there are {count} such groups.")

solve_group_problem()
