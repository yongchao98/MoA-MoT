import itertools

def is_derangement(p):
    """
    Checks if a given 0-indexed permutation is a derangement.
    A permutation is a derangement if p[i] != i for all i.
    """
    for i, pi in enumerate(p):
        if i == pi:
            return False
    return True

def count_excedances(p):
    """
    Counts the number of excedances in a 0-indexed permutation.
    An excedance is an index i such that p[i] > i.
    """
    count = 0
    for i, pi in enumerate(p):
        if pi > i:
            count += 1
    return count

def get_derangement_poly_coeffs(n):
    """
    Calculates the coefficients of the n-th derangement polynomial d_n(t).
    Returns a dictionary where keys are exponents (number of excedances)
    and values are coefficients (number of derangements with that many excedances).
    """
    # For n=0, d_0(t) = 1 (by convention, permutation of empty set)
    # For n=1, d_1(t) = 0 (no derangements)
    if n < 2:
        return {}
        
    coeffs = {}
    elements = range(n)
    for p in itertools.permutations(elements):
        if is_derangement(p):
            exc = count_excedances(p)
            coeffs[exc] = coeffs.get(exc, 0) + 1
    return coeffs

def solve_and_print():
    """
    Solves the three parts of the problem and prints the answer in the specified format.
    """
    # (a) Confirm whether H(U_{n-1, E})(t) = t^{n-1} d_n(t).
    # Based on theoretical degree analysis:
    # deg(H(U_{n-1, E})(t)) = n - 2
    # deg(t^{n-1} d_n(t)) = (n - 1) + deg(d_n(t)) = (n - 1) + (n - 1) = 2n - 2.
    # Since n - 2 != 2n - 2 for n >= 2, the equality is false.
    answer_a = "No"

    # (b) State if the leading coefficient of d_n(t) for any n >= 2 is always 1.
    # The leading coefficient corresponds to the maximum number of excedances, which is n-1.
    # The only permutation with n-1 excedances is the cycle (2, 3, ..., n, 1), which
    # is a derangement. Therefore, its coefficient is 1. This is confirmed for n=2..5.
    is_leading_coeff_one = True
    for n in range(2, 6): # Check for n=2, 3, 4, 5 to be sure
        coeffs = get_derangement_poly_coeffs(n)
        if not coeffs:
            is_leading_coeff_one = False
            break
        max_exc = max(coeffs.keys())
        # The degree of d_n(t) for n>=2 is n-1.
        if max_exc != n - 1 or coeffs[max_exc] != 1:
            is_leading_coeff_one = False
            break
    answer_b = "Yes" if is_leading_coeff_one else "No"

    # (c) Give the value of d_3(1).
    # d_3(1) is the number of derangements of 3 elements.
    n_c = 3
    coeffs_3 = get_derangement_poly_coeffs(n_c)
    # The value of the polynomial at t=1 is the sum of its coefficients.
    num_derangements_3 = sum(coeffs_3.values())
    answer_c = num_derangements_3

    # Print the final combined answer.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_and_print()