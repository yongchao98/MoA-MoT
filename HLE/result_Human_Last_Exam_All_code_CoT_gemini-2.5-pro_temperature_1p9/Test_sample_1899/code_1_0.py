import numpy as np

def find_inf_sup():
    """
    This function calculates the infimum and supremum for the given problem.
    It computes the sequence E_n = |P_n(xi)| * (a_n^2 + b_n^2 + c_n^2) and analyzes its behavior.
    """
    # Define the polynomial f(x) = x^3 + x^2 + x - 1
    p = [1, 1, 1, -1]
    
    # Find the roots of the polynomial
    roots = np.roots(p)
    
    # The real root xi is the first one in the output of numpy.roots
    xi = roots[0].real
    
    # Initial coefficients for n=1 (P_1(x) = x)
    a, b, c = 0, 1, 0
    
    # Store the values of E_n to find min and max
    e_values = []
    
    # Iterate for a sufficient number of steps to observe the trend
    num_iterations = 50 
    
    # We will track the minimum value found and the maximum value
    infimum = float('inf')
    supremum = float('-inf')

    print(f"The real root is xi = {xi:.6f}\n")
    print("Calculating the first few terms of the sequence E_n:")
    
    for n in range(1, num_iterations + 1):
        # P_n(xi) is xi^n. Since xi is positive, |P_n(xi)| = xi^n.
        p_n_xi = xi**n
        
        # Calculate the sum of squares of coefficients
        sum_sq = a**2 + b**2 + c**2
        
        # Calculate E_n
        e_n = p_n_xi * sum_sq
        e_values.append(e_n)

        # Update infimum and supremum
        if e_n < infimum:
            infimum = e_n
        if e_n > supremum:
            supremum = e_n
            
        if n <= 15:
             print(f"For n={n:2d}, (a,b,c) = ({a:4d},{b:4d},{c:4d}), a^2+b^2+c^2 = {sum_sq:8d}, E_{n:2d} = {e_n:.6f}")
        elif n == 16:
             print("...")

        # Update coefficients for the next iteration using the recurrence relation
        # a_{n+1} = c_n, b_{n+1} = a_n - c_n, c_{n+1} = b_n - c_n
        a_next = c
        b_next = a - c
        c_next = b - c
        
        a, b, c = a_next, b_next, c_next

    # Check if the sequence appears to be unbounded
    if supremum > 1e6:
      supremum_val = "infinity"
    else:
      # In case the sequence is bounded, this would be the observed maximum
      supremum_val = supremum
    
    print("\nBased on numerical evidence:")
    print(f"The infimum of |P_n(xi)|(a_n^2+b_n^2+c_n^2) is {infimum}")
    print(f"The supremum of |P_n(xi)|(a_n^2+b_n^2+c_n^2) is {supremum_val}")


find_inf_sup()