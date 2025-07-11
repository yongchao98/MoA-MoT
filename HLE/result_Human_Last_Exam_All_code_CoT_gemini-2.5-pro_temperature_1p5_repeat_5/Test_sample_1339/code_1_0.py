import sys

def solve():
    """
    This function prints the solution to the theoretical group theory problem.
    """
    # Part (a): Based on the theory of localization of solvable groups, a unique minimal
    # group G_hat with the specified properties exists.
    answer_a = "Yes"

    # Part (b): The derived length of G is bounded by the length of its subnormal
    # series with abelian factors, which is n. This bound is achievable, for example,
    # by iterated wreath products of p'-torsion-free abelian groups. The derived
    # length is preserved during p-completion.
    # The variable 'n' is given in the problem statement as the length of the series.
    answer_b = "n"

    # The problem asks for the answer in a specific format.
    print(f"(a) {answer_a}; (b) {answer_b}")

solve()