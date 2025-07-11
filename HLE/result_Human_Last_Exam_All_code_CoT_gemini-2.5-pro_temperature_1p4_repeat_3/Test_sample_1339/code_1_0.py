def solve_group_theory_problem():
    """
    This function provides the solution to the group theory problem.
    The reasoning is as follows:

    Part (a): The group `hat{G}` is the p-localization of G.
    This is a standard construction with a universal property,
    which ensures that the minimal such group is unique up to isomorphism.
    Therefore, the answer is 'Yes'.

    Part (b): The derived length of `hat{G}` is the same as the derived length of G.
    The existence of the given subnormal series of length n with abelian factors
    implies that the derived length of G is at most n.
    This maximum is achievable, for instance, by considering a finite p-group
    of derived length n and its derived series. The factors of this series
    are p-groups and thus p'-torsion-free.
    Therefore, the maximum possible derived length is n.
    """
    answer_a = "Yes"
    answer_b = "n"
    
    print(f"(a) {answer_a}; (b) {answer_b}")

solve_group_theory_problem()