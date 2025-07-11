import sympy
from sympy import symbols, log, Rational, simplify

def calculate_l_k_n():
    """
    Calculates and prints the exact value of l_k(n) in terms of n and k.
    
    The function follows the analytical derivation outlined above:
    l_k(n) = (Jacobian Term) + (Determinant Term) - (Quadratic Form Term)
    """
    # Define the symbolic variables n and k
    n, k = symbols('n k', real=True, positive=True)

    # 1. Term from the Jacobian of the exponential map
    # This comes from (n-1) * ln(||v|| / sinh(||v||))
    # where ||v|| = k and sinh(k) = 1.
    # The term is (n-1) * ln(k)
    term_jacobian = (n - 1) * log(k)
    # The numbers in this part of the equation are:
    num_jac_n_coeff = 1
    num_jac_n_sub = 1
    num_jac_log_coeff = 1
    
    # 2. Term from the determinant of the covariance matrix
    # This comes from -1/2 * ln(det(Sigma))
    # where det(Sigma) = 1/(n+1).
    # The term is -1/2 * ln(1/(n+1)) = 1/2 * ln(n+1)
    term_determinant = Rational(1, 2) * log(n + 1)
    # The numbers in this part of the equation are:
    num_det_coeff_num = 1
    num_det_coeff_den = 2
    num_det_n_coeff = 1
    num_det_n_add = 1

    # 3. Term from the quadratic form in the exponent
    # This comes from 1/2 * n_eval^T * Sigma^-1 * n_eval
    # The result of this calculation is (k^2 * (5*n - 2)) / (2*n)
    term_quadratic = (k**2 * (5*n - 2)) / (2*n)
    # The numbers in this part of the equation are:
    num_quad_k_pow = 2
    num_quad_n_coeff = 5
    num_quad_sub = 2
    num_quad_den_n_coeff = 2
    
    # Combine the terms to get the final expression for l_k(n)
    l_k_n = term_determinant + term_jacobian - term_quadratic

    # Print the final result
    print("The exact value of l_k(n) is:")
    sympy.pprint(l_k_n)
    
    print("\nIn standard notation, the formula is:")
    print("l_k(n) = (1/2)*ln(n+1) + (n-1)*ln(k) - (k**2*(5*n-2))/(2*n)")
    
    # To satisfy the "output each number" requirement, we list the coefficients explicitly.
    print("\nThe numerical constants in the final equation are:")
    print(f"Coefficient of ln(n+1): {num_det_coeff_num}/{num_det_coeff_den}")
    print(f"Term (n - c): c = {num_jac_n_sub}")
    print(f"Coefficient of ln(k): {num_jac_log_coeff}")
    print("In the term -(k**A*(B*n-C))/(D*n):")
    print(f"  Exponent A on k: {num_quad_k_pow}")
    print(f"  Coefficient B of n in numerator: {num_quad_n_coeff}")
    print(f"  Constant C subtracted in numerator: {num_quad_sub}")
    print(f"  Coefficient D of n in denominator: {num_quad_den_n_coeff}")


if __name__ == "__main__":
    calculate_l_k_n()
