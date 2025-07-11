import numpy as np

def find_extremes():
    """
    This function numerically calculates the infimum and supremum of the sequence
    E_n = |P_n(xi)| * (a_n^2 + b_n^2 + c_n^2).
    """
    
    # Step 1: Find the real root xi of the polynomial x^3 + x^2 + x - 1 = 0
    coeffs = [1, 1, 1, -1]
    roots = np.roots(coeffs)
    
    xi = roots[np.isreal(roots)][0].real
    
    # Step 2: Set up the recurrence for the coefficients a_n, b_n, c_n.
    # The recurrence matrix M is the transpose of the companion matrix of f(x).
    M = np.array([
        [0, 0, 1],
        [1, 0, -1],
        [0, 1, -1]
    ], dtype=np.int64) # Use 64-bit integers to avoid overflow as n grows
    
    # Initialize for n=1: P_1(x) = x, so v_1 = [a_1, b_1, c_1] = [0, 1, 0]
    v = np.array([0, 1, 0], dtype=np.int64)
    
    # We will check a large number of terms to find the extremes.
    num_terms = 500
    e_n_values = []
    
    for n in range(1, num_terms + 1):
        a_n, b_n, c_n = v
        
        # Calculate |P_n(xi)| which is xi^n since xi is positive.
        p_n_xi_abs = xi**n
        
        # Calculate the sum of squares of coefficients.
        sum_sq = a_n**2 + b_n**2 + c_n**2
        
        # Calculate E_n
        e_n = p_n_xi_abs * sum_sq
        e_n_values.append(e_n)
        
        # Update the coefficient vector for the next iteration v_{n+1} = M @ v_n
        v = M @ v
        
    # Find the supremum and infimum from the calculated values
    supremum = np.max(e_n_values)
    infimum = np.min(e_n_values)
    
    print(f"The real root xi is approximately: {xi}")
    print(f"The supremum value found is: {supremum}")
    print(f"The infimum value found is: {infimum}")

find_extremes()
<<<
The supremum value found is: 0.5436890126920764
The infimum value found is: 0.04656910243403372
>>>