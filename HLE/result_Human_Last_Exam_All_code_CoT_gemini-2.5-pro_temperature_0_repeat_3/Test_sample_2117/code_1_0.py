def solve_matrix_problem():
    """
    Solves the problem by analyzing the properties of the Gaussian Hessenberg Decomposition.

    The problem asks for the least upper bound of the product E_P * E_H * S_P * S_H.

    1. The decomposition is C_n = P * H * P^-1, where P is derived from
       Gaussian elimination steps.
    2. The standard algorithm produces a unit lower triangular matrix P.
    3. The eigenvalues of a triangular matrix are its diagonal entries.
       For a unit triangular matrix, all diagonal entries are 1.
    4. Thus, all eigenvalues of P are 1.
    5. The average eigenvalue gap E_P is defined as (lambda_max(P) - lambda_min(P)) / (size-1).
       For P, lambda_max(P) = 1 and lambda_min(P) = 1.
    6. This makes E_P = (1 - 1) / (n + 1) = 0 for any positive integer n.
    7. Since E_P is a factor in the product, the entire product is always 0.
    8. The set of values for the product over all n is {0}.
    9. The least upper bound of {0} is 0.
    """

    # For any positive integer n, the matrix P has eigenvalues equal to 1.
    lambda_max_P = 1
    lambda_min_P = 1

    # The size of P is (n+2) x (n+2). The formula for E_P is
    # (lambda_max - lambda_min) / (size - 1).
    # For any n > 0, the denominator (n+1) is positive.
    # The numerator is lambda_max_P - lambda_min_P.
    numerator_E_P = lambda_max_P - lambda_min_P
    
    # E_P is therefore 0 for all n.
    E_P = 0

    print(f"The maximum eigenvalue of P is: {lambda_max_P}")
    print(f"The minimum eigenvalue of P is: {lambda_min_P}")
    print(f"The average eigenvalue gap of P, E_P, is ({lambda_max_P} - {lambda_min_P}) / (n+1) = {E_P}")
    
    # The full product is E_P * E_H * S_P * S_H.
    # Since E_P = 0, the entire product is 0.
    final_product = 0
    
    print(f"The final product is E_P * E_H * S_P * S_H = {E_P} * E_H * S_P * S_H = {final_product}")

    # The product is 0 for all n. The least upper bound of the set {0} is 0.
    least_upper_bound = 0
    print(f"The least upper bound over all positive integers n is: {least_upper_bound}")

solve_matrix_problem()