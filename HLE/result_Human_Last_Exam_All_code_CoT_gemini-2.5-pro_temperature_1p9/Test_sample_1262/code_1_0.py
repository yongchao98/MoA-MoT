import math
from itertools import permutations

def is_derangement(perm_tuple):
    """
    Checks if a permutation is a derangement.
    A permutation p of {0, 1, ..., n-1} is a derangement if p[i] != i for all i.
    """
    for i in range(len(perm_tuple)):
        if perm_tuple[i] == i:
            return False
    return True

def count_excedances(perm_tuple):
    """
    Counts the number of excedances in a permutation.
    An excedance is an index i such that p[i] > i.
    """
    count = 0
    for i in range(len(perm_tuple)):
        if perm_tuple[i] > i:
            count += 1
    return count

def get_derangement_polynomial_coeffs(n):
    """
    Calculates the coefficients of the derangement polynomial d_n(t).
    Returns a dictionary mapping exponent to coefficient.
    """
    coeffs = {}
    if n == 0:
        return {0: 1} # by convention
    if n == 1:
        return {} # d_1(t) = 0
        
    # Generate permutations of {0, 1, ..., n-1}
    all_perms = permutations(range(n))
    
    for p in all_perms:
        if is_derangement(p):
            exc_count = count_excedances(p)
            coeffs[exc_count] = coeffs.get(exc_count, 0) + 1
    return coeffs

def solve():
    """
    Solves all parts of the user's question.
    """
    # --- Part (a) ---
    # Test the identity H(U_{n-1, E})(t) = t^(n-1) * d_n(t) for n=3.
    n_test = 3
    # LHS: Hilbert series of the Chow ring of U_{n-1, n}
    # For n>=2, H(t) = 1 + (C(n, n-2) - 1)*t
    h_poly_coeff_t1 = math.comb(n_test, n_test - 2) - 1
    lhs_poly_coeffs = {0: 1, 1: h_poly_coeff_t1}

    # RHS: t^(n-1) * d_n(t)
    d_n_coeffs = get_derangement_polynomial_coeffs(n_test)
    rhs_poly_coeffs = {k + (n_test - 1): v for k, v in d_n_coeffs.items()}
    
    # Compare LHS and RHS
    # For n=3, LHS is 1 + 2t. RHS is t^2 * (t + t^2) = t^3 + t^4. They are not equal.
    ans_a = "No"

    # --- Part (b) ---
    # Check the leading coefficient of d_n(t) for n >= 2.
    # The leading coefficient is the number of derangements with max excedances.
    # The max number of excedances is n-1, achieved only by sigma = (2,3,...,n,1).
    # Thus, the leading coefficient is 1. We verify for a few n.
    is_b_true = True
    for n in range(2, 6):
        d_coeffs = get_derangement_polynomial_coeffs(n)
        if d_coeffs:
            max_power = max(d_coeffs.keys())
            leading_coeff = d_coeffs[max_power]
            if max_power != n - 1 or leading_coeff != 1:
                is_b_true = False
                break
    
    if is_b_true:
        ans_b = "Yes"
    else:
        ans_b = "No"
    
    # --- Part (c) ---
    # Calculate d_3(1), which is the number of derangements of 3 elements.
    d3_coeffs = get_derangement_polynomial_coeffs(3)
    # d_3(1) is the sum of coefficients.
    ans_c = sum(d3_coeffs.values())

    # --- Format final answer ---
    print(f"(a) {ans_a}; (b) {ans_b}; (c) {ans_c}")

solve()