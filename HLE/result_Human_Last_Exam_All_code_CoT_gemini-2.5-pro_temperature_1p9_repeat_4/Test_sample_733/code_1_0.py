def solve_group_theory_question():
    """
    This function provides the number of finite groups containing
    a maximal by inclusion product-free set of size 2.

    This solution is based on a known classification theorem from combinatorial group theory.
    A subset S of a group G is product-free if for all x, y in S, the product xy is not in S.
    A product-free set S is maximal if for any element g in G \ S, the set S U {g} is not product-free.

    The problem asks for the number of non-isomorphic finite groups that have such a set S of size 2.
    According to the paper "Small maximal product-free sets" by Giudici and Hart (2009),
    there are exactly 8 such groups.
    """

    # The list of non-isomorphic finite groups with a maximal product-free set of size 2
    groups = [
        "C_3 (Cyclic group of order 3)",
        "C_4 (Cyclic group of order 4)",
        "C_2 x C_2 (Klein four-group)",
        "C_5 (Cyclic group of order 5)",
        "S_3 (Symmetric group of order 6, also known as D_3)",
        "C_7 (Cyclic group of order 7)",
        "A_4 (Alternating group of order 12)",
        "Q_12 (Dicyclic group of order 12)"
    ]

    print("The finite groups containing a maximal by inclusion product-free set of size 2 are:")
    for group in groups:
        print(f"- {group}")

    count = len(groups)
    print(f"\nThere are a total of {count} such groups.")

if __name__ == "__main__":
    solve_group_theory_question()