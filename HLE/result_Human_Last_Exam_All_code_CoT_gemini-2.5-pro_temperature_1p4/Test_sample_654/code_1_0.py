import sympy

def solve_polynomial_problem(p, k, n):
    """
    This function computes the number of coefficients in the given recursively defined polynomial
    that are not divisible by p^k.

    Args:
      p: An odd prime.
      k: A positive integer.
      n: A positive integer.

    Returns:
      The number of coefficients not divisible by p^k.
    """
    
    # The reasoning for the solution is as follows:
    # Let P_j(x) be the polynomial after j iterations. We want to find the number of
    # coefficients of P_{p^n}(x) that are not divisible by p^k.
    # Let's analyze the polynomials modulo p. Let Q_j(x) = P_j(x) mod p.
    # The recurrence relation is Q_{j+1}(x) = (Q_j(x)^p - Q_j(x))^(p^(k-1)) mod p.
    #
    # Q_0(x) = x
    # Q_1(x) = (x^p - x)^(p^(k-1)) = x^(p^k) - x^(p^(k-1)) mod p. (2 terms)
    #
    # Let R_j(x) = Q_{j-1}(x)^p - Q_{j-1}(x). The number of non-zero coefficients in Q_j is the same
    # as in R_j.
    # The set of exponents in R_j(x) forms a clean pattern.
    # Let's denote the number of non-zero coefficients (terms) in Q_j(x) as T(Q_j).
    #
    # T(Q_0) = 1.
    # R_1(x) = Q_0^p - Q_0 = x^p - x. Has 2 terms.
    # T(Q_1) = 2.
    #
    # R_2(x) = Q_1^p - Q_1 = (x^(p^k) - x^(p^(k-1)))^p - (x^(p^k) - x^(p^(k-1)))
    #        = (x^(p^(k+1)) - x^(p^k)) - (x^(p^k) - x^(p^(k-1)))
    #        = x^(p^(k+1)) - 2x^(p^k) + x^(p^(k-1)) mod p. Has 3 terms.
    # T(Q_2) = 3.
    #
    # R_3(x) = Q_2^p - Q_2. The exponents in Q_2 are p^(2k-2), p^(2k-1), p^(2k).
    # The exponents in Q_2^p are p^(2k-1), p^(2k), p^(2k+1).
    # The union of these two sets has size 4.
    # T(Q_3) = 4.
    #
    # By induction, the number of non-zero coefficients of Q_j(x) mod p is j+1.
    # For the j-th step being p^n, the number of non-zero coefficients is p^n + 1.
    # These coefficients are not divisible by p, so they are not divisible by p^k.
    # This argument suggests that the answer is p^n + 1.
    
    result = p**n + 1
    
    # Printing the result in the equation form
    print(f"p = {p}")
    print(f"k = {k}")
    print(f"n = {n}")
    print(f"The number of coefficients not divisible by p^k is given by the expression: p^n + 1")
    print(f"Substituting the values: {p}^{n} + 1 = {p**n} + 1 = {result}")
    
    return result

# Example usage (not part of the final output, just for testing)
if __name__ == '__main__':
    # This part is for testing the logic, but the user expects the function block.
    # Let's use some example values. p=3, k=1, n=1
    p_val = 3
    k_val = 2
    n_val = 2
    solve_polynomial_problem(p_val, k_val, n_val)