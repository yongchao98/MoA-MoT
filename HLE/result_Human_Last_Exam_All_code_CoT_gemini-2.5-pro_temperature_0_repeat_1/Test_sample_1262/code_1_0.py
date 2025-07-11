import itertools

def solve_derangement_problem():
    """
    Solves the three-part problem concerning derangement polynomials and matroid theory.
    """

    # Part (a): Confirm whether H(U_{n-1, E})(t) = t^{n-1} d_n(t).
    #
    # The Hilbert series H(M)(t) of a matroid M of rank r is a polynomial of degree r.
    # For the uniform matroid U_{n-1, E}, the rank is r = n-1.
    # Therefore, the degree of H(U_{n-1, E})(t) is n-1.
    #
    # The derangement polynomial d_n(t) is the sum of t^exc(sigma) over all derangements.
    # The maximum number of excedances in a derangement of n elements is n-1 (achieved by the
    # permutation (2, 3, ..., n, 1)). Thus, the degree of d_n(t) is n-1 for n >= 2.
    #
    # The degree of the right-hand side, t^{n-1} * d_n(t), would be:
    # deg(t^{n-1}) + deg(d_n(t)) = (n-1) + (n-1) = 2n - 2.
    #
    # For the equality to hold, the degrees must be equal:
    # n - 1 = 2n - 2  => n = 1.
    # The identity fails for any n >= 2. The established identity is actually H(U_{n-1, E})(t) = d_n(t).
    # Therefore, the statement is false.
    answer_a = "No"

    # Part (b): State if the leading coefficient of d_n(t) for any n >= 2 is always 1.
    #
    # The leading coefficient is the number of derangements with the maximum number of excedances.
    # As established, the maximum number of excedances is n-1.
    # We need to find how many derangements sigma have exc(sigma) = n-1.
    # A permutation with n-1 excedances must satisfy sigma(i) > i for n-1 indices i.
    # The only index k where this can fail is k=n, since sigma(n) cannot be > n.
    # So, for a derangement sigma with n-1 excedances, we must have:
    # 1. sigma(i) > i for i = 1, ..., n-1
    # 2. sigma(n) < n (since it's a derangement and not an excedance)
    # From (1), we have sigma(n-1) > n-1, which implies sigma(n-1) = n.
    # Then sigma(n-2) > n-2 and sigma(n-2) != n, so sigma(n-2) must be n-1.
    # Continuing this logic, we find sigma(i) = i+1 for i = 1, ..., n-1.
    # This forces the only remaining value for sigma(n) to be 1.
    # The unique permutation is (2, 3, ..., n, 1). This is a derangement for n >= 2.
    # Since there is only one such derangement, the leading coefficient is 1.
    answer_b = "Yes"

    # Part (c): Give the value of d_3(1).
    #
    # d_n(1) is the value of the polynomial when t=1.
    # d_n(1) = sum_{sigma is a derangement} 1^exc(sigma) = number of derangements of n elements.
    # We need to find the number of derangements in S_3.
    # The permutations of {1, 2, 3} are:
    # (1, 2, 3) - not a derangement (3 fixed points)
    # (1, 3, 2) - not a derangement (1 is a fixed point)
    # (2, 1, 3) - not a derangement (3 is a fixed point)
    # (2, 3, 1) - is a derangement
    # (3, 1, 2) - is a derangement
    # (3, 2, 1) - not a derangement (2 is a fixed point)
    # There are 2 derangements.
    
    n_c = 3
    perms_n_c = list(itertools.permutations(range(1, n_c + 1)))
    
    def is_derangement(p):
        for i, val in enumerate(p):
            if i + 1 == val:
                return False
        return True

    num_derangements_c = sum(1 for p in perms_n_c if is_derangement(p))
    answer_c = num_derangements_c

    # Format and print the final answer.
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer)

solve_derangement_problem()