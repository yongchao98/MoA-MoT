import numpy as np
from decimal import Decimal, getcontext

def solve_numerical_problem():
    """
    This function numerically computes the infimum and supremum of the expression
    |P_n(xi)| * (a_n^2 + b_n^2 + c_n^2) over positive integers n.

    The plan is as follows:
    1. Find the real root xi of the polynomial f(x) = x^3 + x^2 + x - 1 with high precision.
    2. Use the matrix recurrence relation for the coefficients (a_n, b_n, c_n) to generate them iteratively.
    3. For each n up to a large limit, calculate the expression E_n = xi^n * (a_n^2 + b_n^2 + c_n^2).
    4. Keep track of the minimum and maximum values of E_n encountered.
    5. Print the final results, including the details of the calculation for the n where the infimum and supremum were found.
    """
    # Set precision for Decimal calculations
    getcontext().prec = 50

    # Step 1: Find the real root xi with high precision using Newton's method
    # Use numpy to get a good initial guess
    initial_guess = np.roots([1, 1, 1, -1])
    real_root_guess = [r.real for r in initial_guess if abs(r.imag) < 1e-9][0]
    
    x = Decimal(real_root_guess)
    f = lambda t: t**3 + t**2 + t - 1
    f_prime = lambda t: 3*t**2 + 2*t + 1
    
    # Perform a few iterations of Newton's method to refine the root
    for _ in range(10):
        x = x - f(x) / f_prime(x)
    xi = x
    
    # Step 2: Set up the recurrence
    # Initial vector for n=1: xi^1 = 0 + 1*xi + 0*xi^2 => (a_1,b_1,c_1) = (0,1,0)
    V_n = np.array([0, 1, 0], dtype=np.longlong)
    
    # Recurrence matrix M
    M = np.array([[0, 0, 1], 
                  [1, 0, -1], 
                  [0, 1, -1]], dtype=np.longlong)

    min_E = Decimal('inf')
    max_E = Decimal('-inf')
    
    # Set a limit for n. A few thousand should be sufficient to find the extrema.
    n_max = 4000
    
    # Pre-calculate powers of xi to avoid repeated multiplications
    xi_powers = [xi]
    for _ in range(n_max - 1):
        xi_powers.append(xi_powers[-1] * xi)

    # Variables to store details for the final printout
    inf_n, sup_n = 0, 0
    inf_coeffs, sup_coeffs = None, None
    inf_Pn_xi, sup_Pn_xi = None, None

    # Step 3: Iterate and calculate E_n
    for n in range(1, n_max + 1):
        if n > 1:
            V_n = M @ V_n
            
        a_n, b_n, c_n = V_n[0], V_n[1], V_n[2]
        
        # Calculate sum of squares
        sum_sq = Decimal(a_n**2 + b_n**2 + c_n**2)
        
        # |P_n(xi)| is xi^n since xi is positive
        Pn_xi = xi_powers[n-1]
        
        # The expression to evaluate
        E_n = Pn_xi * sum_sq
        
        # Step 4: Keep track of min and max values and their details
        if E_n < min_E:
            min_E = E_n
            inf_n = n
            inf_coeffs = (a_n, b_n, c_n)
            inf_Pn_xi = Pn_xi

        if E_n > max_E:
            max_E = E_n
            sup_n = n
            sup_coeffs = (a_n, b_n, c_n)
            sup_Pn_xi = Pn_xi
            
    # Step 5: Output the results
    print("# Finding the infimum")
    a_inf, b_inf, c_inf = inf_coeffs
    sum_sq_inf = a_inf**2 + b_inf**2 + c_inf**2
    print(f"The infimum is approximately {min_E:.10f}, found at n = {inf_n}")
    print(f"The calculation for the infimum is:")
    print(f"|P_{inf_n}(xi)| * (a_{inf_n}^2 + b_{inf_n}^2 + c_{inf_n}^2) = {inf_Pn_xi:.10f} * ({a_inf}^2 + {b_inf}^2 + {c_inf}^2)")
    print(f"= {inf_Pn_xi:.10f} * {sum_sq_inf} = {min_E:.10f}")
    
    print("\n# Finding the supremum")
    a_sup, b_sup, c_sup = sup_coeffs
    sum_sq_sup = a_sup**2 + b_sup**2 + c_sup**2
    print(f"The supremum is approximately {max_E:.10f}, found at n = {sup_n}")
    print(f"The calculation for the supremum is:")
    print(f"|P_{sup_n}(xi)| * (a_{sup_n}^2 + b_{sup_n}^2 + c_{sup_n}^2) = {sup_Pn_xi:.10f} * ({a_sup}^2 + {b_sup}^2 + {c_sup}^2)")
    print(f"= {sup_Pn_xi:.10f} * {sum_sq_sup} = {max_E:.10f}")

if __name__ == '__main__':
    solve_numerical_problem()