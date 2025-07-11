import numpy as np

def solve():
    """
    This function calculates the infimum and supremum of the expression
    |P_n(xi)| * (a_n^2 + b_n^2 + c_n^2) for the given polynomial.
    """
    # 1. Define the polynomial and find its real root xi.
    # The polynomial is f(x) = x^3 + x^2 + x - 1
    coeffs = [1, 1, 1, -1]
    roots = np.roots(coeffs)
    
    # Find the real root xi
    xi = None
    for r in roots:
        if np.isreal(r):
            xi = r.real
            break
    
    # 2. Initialize coefficients for n=1
    # For n=1, P_1(x) = x, so a_1=0, b_1=1, c_1=0.
    a, b, c = 0, 1, 0

    # 3. Calculate the first value of the expression for n=1
    # E_n = xi^n * (a_n^2 + b_n^2 + c_n^2)
    norm_sq = a**2 + b**2 + c**2
    e_n = (xi**1) * norm_sq
    
    infimum = e_n
    supremum = e_n

    # 4. Iterate for n from 2 to a large number (e.g., 100)
    num_iterations = 100
    for n in range(2, num_iterations + 1):
        # Update coefficients using the recurrence relation:
        # a_{n+1} = c_n
        # b_{n+1} = a_n - c_n
        # c_{n+1} = b_n - c_n
        a_prev, b_prev, c_prev = a, b, c
        a = c_prev
        b = a_prev - c_prev
        c = b_prev - c_prev
        
        # Calculate the expression E_n
        norm_sq = a**2 + b**2 + c**2
        e_n = (xi**n) * norm_sq
        
        # Update infimum and supremum
        if e_n < infimum:
            infimum = e_n
        if e_n > supremum:
            supremum = e_n

    print(f"The real root xi is: {xi}")
    print(f"Infimum found: {infimum}")
    print(f"Supremum found: {supremum}")

solve()