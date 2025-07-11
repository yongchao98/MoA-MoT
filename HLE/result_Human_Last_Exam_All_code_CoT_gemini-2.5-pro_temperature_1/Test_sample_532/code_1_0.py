def solve_finite_filled_nilpotent_groups():
    """
    This function provides the classification of finite filled nilpotent groups.

    A finite group G is called a filled group if the union of all its maximal
    by inclusion product-free sets is equal to G \ {e}.

    A finite group is nilpotent if and only if it is the direct product of its
    Sylow p-subgroups.
    """

    answer_line_1 = "A finite nilpotent group G is a direct product of its Sylow p-subgroups."
    answer_line_2 = "G is filled if and only if each of its Sylow p-subgroups is filled."
    answer_line_3 = "The classification of filled p-groups is as follows:"
    answer_line_4 = "1. If p is an odd prime, any finite p-group is filled."
    answer_line_5 = "2. If p = 2, a finite 2-group P is filled if and only if its exponent is at most 4."
    answer_line_6 = "(The exponent of a group is the smallest integer n > 0 such that g^n = e for all g in P)."

    conclusion_header = "\nConclusion:"
    conclusion = (
        "A finite nilpotent group G is a filled group if and only if its "
        "Sylow 2-subgroup P_2 has an exponent of at most 4."
    )
    conclusion_detail = "If G has an odd order (i.e., no Sylow 2-subgroup), the condition is trivially met, and the group is filled."

    print(answer_line_1)
    print(answer_line_2)
    print(answer_line_3)
    print("    " + answer_line_4)
    print("    " + answer_line_5)
    print("    " + answer_line_6)
    print(conclusion_header)
    print(conclusion)
    print(conclusion_detail)

solve_finite_filled_nilpotent_groups()