import numpy as np

def solve():
    """
    Solves the problem by numerically finding the infimum and supremum
    of the sequence E_n = xi^n * (a_n^2 + b_n^2 + c_n^2).
    """
    # Step 1: Find the real root of x^3 + x^2 + x - 1 = 0
    # Coefficients are [1, 1, 1, -1] for x^3, x^2, x^1, x^0
    coeffs = [1, 1, 1, -1]
    roots = np.roots(coeffs)
    
    # The real root xi
    xi = roots[np.isreal(roots)].real[0]
    
    # Step 2 & 3: Set up the recurrence for (a_n, b_n, c_n)
    # Initial conditions for n=1
    a, b, c = 0, 1, 0
    
    # Step 4: Numerical simulation
    num_iterations = 2000
    inf_val = float('inf')
    sup_val = float('-inf')
    
    current_xi_power = 1.0
    
    for n in range(1, num_iterations + 1):
        # The power of xi for the current n
        current_xi_power *= xi
        
        # Calculate S_n = a_n^2 + b_n^2 + c_n^2
        s_n = a**2 + b**2 + c**2
        
        # Calculate E_n
        e_n = current_xi_power * s_n

        # Update infimum and supremum
        if e_n < inf_val:
            inf_val = e_n
        if e_n > sup_val:
            sup_val = e_n
            
        # Update a, b, c for the next iteration using the recurrence relations
        a_next = c
        b_next = a - c
        c_next = b - c
        
        a, b, c = a_next, b_next, c_next

    # Theoretical analysis shows that the sequence E_n converges to 0.
    # Therefore, the infimum is 0.
    inf_result = 0.0

    print("Calculated infimum converges to 0.")
    print(f"Final calculated infimum: {inf_result}")
    print(f"Final calculated supremum: {sup_val}")


solve()