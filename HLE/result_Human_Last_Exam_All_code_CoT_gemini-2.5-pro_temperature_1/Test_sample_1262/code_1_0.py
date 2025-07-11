import itertools

def get_derangement_polynomial_coeffs(n):
    """
    Computes the coefficients of the n-th derangement polynomial d_n(t).
    Returns a list of coefficients, where the index corresponds to the power of t.
    e.g., [c_0, c_1, c_2] for c_0 + c_1*t + c_2*t^2
    """
    if n == 0:
        # By convention, D_0 = 1, so d_0(t) = 1
        return [1]
    if n == 1:
        # No derangements in S_1, so d_1(t) = 0
        return []

    # Maximum excedances in a derangement is n-1
    max_exc = n - 1
    coeffs = [0] * (max_exc + 1)

    # Iterate through all permutations of {1, 2, ..., n}
    for p in itertools.permutations(range(1, n + 1)):
        # Check if it's a derangement
        is_derangement = True
        for i in range(n):
            if p[i] == i + 1:  # Check for fixed points sigma(i) = i
                is_derangement = False
                break
        
        if is_derangement:
            # Count excedances: sigma(i) > i
            excedance_count = 0
            for i in range(n):
                if p[i] > i + 1:
                    excedance_count += 1
            
            if excedance_count <= max_exc:
                coeffs[excedance_count] += 1
    
    return coeffs

def solve():
    """
    Solves the three parts of the problem and prints the formatted answer.
    """
    # Part (a): The identity is H(U_{n-1,E})(t) = t^{n-1}d_n(t).
    # The degree of the LHS is n-1.
    # The degree of the RHS is (n-1) + deg(d_n(t)).
    # For n>=2, deg(d_n(t)) = n-1. So deg(RHS) = 2n-2.
    # n-1 = 2n-2 implies n=1. The identity is false for n>=2.
    answer_a = "No"

    # Part (b): Leading coefficient of d_n(t) for n>=2.
    is_leading_coeff_one = True
    for n in range(2, 6): # Check for a few values of n
        coeffs = get_derangement_polynomial_coeffs(n)
        # Find the last non-zero coefficient
        leading_coeff = 0
        if coeffs:
            for c in reversed(coeffs):
                if c != 0:
                    leading_coeff = c
                    break
        if leading_coeff != 1:
            is_leading_coeff_one = False
            break
    answer_b = "Yes" if is_leading_coeff_one else "No"

    # Part (c): Value of d_3(1).
    # d_3(1) is the number of derangements of 3 elements.
    d3_coeffs = get_derangement_polynomial_coeffs(3)
    num_derangements_3 = sum(d3_coeffs)
    answer_c = num_derangements_3

    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve()