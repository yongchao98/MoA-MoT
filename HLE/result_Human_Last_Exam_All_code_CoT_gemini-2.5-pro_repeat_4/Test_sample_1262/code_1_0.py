import itertools

def is_derangement(p):
    """
    Checks if a 0-indexed permutation tuple 'p' is a derangement.
    A permutation is a derangement if p[i] != i for all i.
    """
    for i, val in enumerate(p):
        if i == val:
            return False
    return True

def count_excedances(p):
    """
    Counts the number of excedances in a 0-indexed permutation tuple 'p'.
    An excedance at index i is when p[i] > i.
    """
    count = 0
    for i, val in enumerate(p):
        if val > i:
            count += 1
    return count

def get_dn_coeffs(n):
    """
    Calculates the coefficients of the n-th derangement polynomial d_n(t).
    The coefficient of t^k is the number of derangements with k excedances.
    Returns a list of coefficients [c_0, c_1, ..., c_{n-1}].
    """
    if n < 0:
        return []
    # The degree of d_n(t) is at most n-1.
    coeffs = [0] * n
    
    # Generate all permutations of {0, 1, ..., n-1}
    for p in itertools.permutations(range(n)):
        if is_derangement(p):
            exc = count_excedances(p)
            if exc < n:
                coeffs[exc] += 1
    return coeffs

def solve_problem():
    """
    Solves all parts of the problem and prints the final formatted answer.
    """
    # Part (a): Based on degree analysis, the identity is false.
    answer_a = "No"

    # Part (b): Check if the leading coefficient of d_n(t) is 1 for n>=2.
    # The degree is n-1. The leading coefficient is coeffs[n-1].
    is_leading_coeff_one = True
    for n in range(2, 8):  # Check for a representative range of n
        coeffs = get_dn_coeffs(n)
        # For n>=2, the degree is n-1, so the leading coeff is the last element.
        if len(coeffs) > 0 and coeffs[n - 1] != 1:
            is_leading_coeff_one = False
            break
    answer_b = "Yes" if is_leading_coeff_one else "No"

    # Part (c): Calculate d_3(1).
    # d_n(1) is the sum of coefficients, which is the number of derangements.
    coeffs_3 = get_dn_coeffs(3)
    answer_c = sum(coeffs_3)

    # Print the final answer in the specified format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_problem()