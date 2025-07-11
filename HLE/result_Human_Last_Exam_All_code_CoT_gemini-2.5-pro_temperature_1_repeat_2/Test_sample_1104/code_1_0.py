def solve_proportionality_problem():
    """
    This function calculates the smallest preference profile sizes s1 and s2.

    s1: Smallest size for Proportional Justified Representation (PJR)
    s2: Smallest size for Extended Justified Representation (EJR)

    The problem is solved analytically as described in the text.
    The committee size k is 100.
    """

    k = 100

    # For both PJR and EJR, a lower bound on the number of voters n can be
    # found by considering the group N' = {voter 1}.
    # Voter 1 must be unsatisfied, meaning A(1) intersect W is empty.
    # PJR/EJR state that if N' is cohesive enough, A(1) intersect W cannot be empty.
    # To avoid this contradiction, N' must NOT be cohesive enough.
    # This leads to the inequality: |N'| < n/k
    # 1 < n/100  => n > 100.
    # The smallest integer n satisfying this is 101.
    s1_lower_bound = k + 1
    s2_lower_bound = k + 1

    # We then show that n=101 is achievable for both cases by constructing a valid
    # profile extension and a winning committee W. The text above provides this
    # constructive proof.
    s1 = s1_lower_bound
    s2 = s2_lower_bound

    # The result is the pair (s1, s2)
    result = (s1, s2)
    print(f"s1 = {result[0]}")
    print(f"s2 = {result[1]}")
    print(f"The solution is the pair (s1, s2): {result}")

solve_proportionality_problem()