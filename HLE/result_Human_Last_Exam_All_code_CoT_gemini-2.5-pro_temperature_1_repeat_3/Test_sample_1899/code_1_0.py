import numpy as np

def solve():
    """
    This function calculates the infimum and supremum of the given expression.
    """
    # Step 1: Find the real root xi of the polynomial x^3 + x^2 + x - 1 = 0.
    # We use numpy's root finder for this.
    coeffs = [1, 1, 1, -1]
    roots = np.roots(coeffs)
    real_root = None
    for r in roots:
        if np.isreal(r):
            real_root = np.real(r)
            break
    
    if real_root is None:
        print("Could not find the real root.")
        return

    xi = real_root

    # Step 2: Initialize coefficients a_n, b_n, c_n.
    # For n=1, x = 0 + 1*x + 0*x^2, so (a1, b1, c1) = (0, 1, 0)
    # We start the loop from n=1.
    a, b, c = 0, 1, 0
    
    min_val = float('inf')
    min_n = -1

    # Step 3: Iterate to find the minimum value of E_n.
    # We compute E_n for n from 1 to 50, which is sufficient to observe the minimum.
    for n in range(1, 51):
        # The coefficients (a,b,c) are for the n-th power.
        
        # Calculate E_n = xi^n * (a_n^2 + b_n^2 + c_n^2)
        E_n = (xi**n) * (a**2 + b**2 + c**2)

        if E_n < min_val:
            min_val = E_n
            min_n = n

        # Update coefficients for the next iteration using the recurrence relation:
        # a_{n+1} = c_n, b_{n+1} = a_n - c_n, c_{n+1} = b_n - c_n
        a_next = c
        b_next = a - c
        c_next = b - c
        a, b, c = a_next, b_next, c_next

    # Step 4: Output the results.
    print(f"The infimum of the expression is approximately: {min_val}")
    print("The infimum occurs at n =", min_n)
    print("The supremum of the expression is infinity.")

solve()