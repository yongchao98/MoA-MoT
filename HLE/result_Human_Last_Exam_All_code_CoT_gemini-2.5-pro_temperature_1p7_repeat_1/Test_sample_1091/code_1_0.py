import sympy
from sympy import Matrix, exp, pi, oo, Symbol, sqrt, pretty_print, Limit, Eq

def calculate_limit():
    """
    This function carries out the plan to find the limit of n*P(n).
    """

    # --- Setup and Definitions ---
    print("This script calculates the limit of n*P(n) as n approaches infinity.")
    print("Plan:")
    print("1. Define the vectors and the number of vectors n.")
    print("2. Calculate the covariance matrix of the sum S using the Central Limit Theorem.")
    print("3. Derive the expression for the probability P(n) by integrating the resulting Gaussian PDF.")
    print("4. Compute the final limit of n*P(n).\n")

    # Symbols for our calculations
    n = Symbol('n', positive=True, real=True)
    r = Symbol('r', positive=True, real=True)
    theta = Symbol('theta', real=True)

    # --- Step 1 & 2: Vectors and Covariance Matrix ---
    print("--- Step 1 & 2: Covariance Matrix Calculation ---")
    
    # Define the three types of vectors
    v_A = Matrix([1, 0])
    v_B = Matrix([sympy.S(1)/2, sqrt(3)/2])
    v_C = Matrix([-sympy.S(1)/2, sqrt(3)/2])

    # The sum of outer products v*v.T for one of each vector
    sum_vvT = v_A * v_A.T + v_B * v_B.T + v_C * v_C.T
    
    # The total number of vectors is n, with n/3 of each type.
    # The covariance matrix of the sum S is Sigma_n = (n/3) * sum_vvT.
    Sigma_n = (n/3) * sum_vvT
    
    print("The covariance matrix of the sum S, denoted as Sigma_n, is:")
    pretty_print(Eq(Symbol('Sigma_n'), Sigma_n))
    print("\n")

    # --- Step 3: Probability P(n) ---
    print("--- Step 3: Probability P(n) Calculation ---")
    
    # For large n, S is approximately N(0, Sigma_n).
    # The PDF of S=(x,y) is f(x,y) = 1 / (2*pi*sqrt(det(Sigma_n))) * exp(-1/2 * [x,y]^T * Sigma_n_inv * [x,y])
    # With Sigma_n = [[n/2, 0], [0, n/2]], the PDF is f(x,y) = 1/(pi*n) * exp(-(x^2+y^2)/n).
    
    # P(n) = Integral of f(x,y) over the disk x^2+y^2 <= 2.
    # In polar coordinates (r, theta), this is Integral[0 to 2*pi] Integral[0 to sqrt(2)] (PDF * r) dr dtheta
    
    # PDF in polar coordinates
    pdf_polar = (1 / (pi * n)) * exp(-r**2 / n)
    
    # We perform the integration to find P(n)
    # The integral with respect to r (including the Jacobian r)
    integral_r = sympy.integrate(r * pdf_polar, (r, 0, sqrt(2)))
    # The integral with respect to theta
    integral_theta = sympy.integrate(integral_r, (theta, 0, 2*pi))
    
    P_n_expr = integral_theta
    
    print("The probability P(n) for large n is approximated by the integral of the Gaussian PDF.")
    print("The resulting expression for P(n) is:")
    pretty_print(Eq(Symbol('P(n)'), P_n_expr))
    print("\n")

    # --- Step 4: Final Limit ---
    print("--- Step 4: Final Limit Calculation ---")
    
    # We want to compute the limit of n * P(n) as n -> oo.
    limit_expr = n * P_n_expr
    limit_value = sympy.limit(limit_expr, n, oo)

    # Display the final equation and its solution
    final_limit_equation = Eq(Limit(limit_expr, n, oo), limit_value)
    
    print("The limit equation is:")
    pretty_print(final_limit_equation)
    
    print("\n----------------------------------------------------")
    print("The final result for the limit of n*P(n) is:")
    print(limit_value)
    print("----------------------------------------------------")

calculate_limit()