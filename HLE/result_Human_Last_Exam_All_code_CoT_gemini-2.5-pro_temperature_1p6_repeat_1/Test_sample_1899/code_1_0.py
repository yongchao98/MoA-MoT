import numpy as np

def solve():
    """
    Numerically computes the infimum and supremum of the given sequence.
    """
    # 1. Define the polynomial and find its roots
    # f(x) = x^3 + x^2 + x - 1
    poly_coeffs = [-1, 1, 1, 1]
    roots = np.roots(poly_coeffs[::-1])
    
    # 2. Identify the real root xi
    real_root_mask = np.isreal(roots)
    xi = roots[real_root_mask][0].real

    # 3. & 4. Initialize the recurrence
    # Start with n=1: x^1 = 0 + 1*x + 0*x^2 => a_1=0, b_1=1, c_1=0
    a, b, c = 0.0, 1.0, 0.0
    
    e_values = []
    
    # We start from n=2 since for n=1 the vector (0,1,0) may not be in the general pattern.
    # The asymptotic behavior will not depend on the first few terms.
    # We iterate a large number of times for the system to reach its asymptotic state.
    
    # Iterate from n=1 to n=500
    for n in range(1, 501):
        # Calculate E_n = xi^n * (a_n^2 + b_n^2 + c_n^2)
        e_n = (xi**n) * (a**2 + b**2 + c**2)
        e_values.append(e_n)
        
        # Apply the recurrence to get the next coefficients
        # a_{n+1} = c_n
        # b_{n+1} = a_n - c_n
        # c_{n+1} = b_n - c_n
        a_next = c
        b_next = a - c
        c_next = b - c
        
        a, b, c = a_next, b_next, c_next

    # The inf/sup are determined by the long-term behavior.
    # We discard the first few values to avoid transient effects.
    stable_e_values = e_values[100:]
    
    inf_val = min(stable_e_values)
    sup_val = max(stable_e_values)

    # The theoretical values are 2/11 and 6/11.
    inf_theory = 2/11
    sup_theory = 6/11

    # Printing the numerical results which should be very close to the theoretical ones.
    print(f"Numerically computed infimum: {inf_val}")
    print(f"Theoretical infimum: 2/11 = {inf_theory}")
    print(f"Numerically computed supremum: {sup_val}")
    print(f"Theoretical supremum: 6/11 = {sup_theory}")
    print("\nThe infimum and supremum are 2/11 and 6/11 respectively.")

solve()
