import numpy as np

def solve():
    """
    This function calculates the infimum and supremum of the expression
    |P_n(xi)| * (a_n^2 + b_n^2 + c_n^2) over positive integers n.
    """
    # 1. Find the real root xi of the polynomial x^3 + x^2 + x - 1 = 0
    coeffs = [1, 1, 1, -1]
    roots = np.roots(coeffs)
    real_root = None
    for root in roots:
        if np.isreal(root):
            real_root = np.real(root)
            break
    
    if real_root is None:
        print("Could not find the real root.")
        return

    xi = real_root

    # 2. Set up initial conditions for n=1
    # For n=1, xi^1 = a_1 + b_1*xi + c_1*xi^2, so a_1=0, b_1=1, c_1=0
    a, b, c = 0.0, 1.0, 0.0

    inf_val = float('inf')
    sup_val = float('-inf')

    # 3. Iterate from n=1 to a sufficiently large number (e.g., 200)
    num_iterations = 200
    for n in range(1, num_iterations + 1):
        # a. Calculate the expression E_n
        xi_n = xi**n
        norm_sq = a**2 + b**2 + c**2
        e_n = xi_n * norm_sq

        # b. Update infimum and supremum
        if e_n < inf_val:
            inf_val = e_n
        if e_n > sup_val:
            sup_val = e_n

        # c. Update coefficients for the next iteration using the recurrence relation
        # a_{n+1} = c_n
        # b_{n+1} = a_n - c_n
        # c_{n+1} = b_n - c_n
        a_next = c
        b_next = a - c
        c_next = b - c
        
        a, b, c = a_next, b_next, c_next

    print(f"Infimum: {inf_val}")
    print(f"Supremum: {sup_val}")

solve()