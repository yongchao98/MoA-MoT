import sympy as sp
import numpy as np
from sympy.functions.special.elliptic_integrals import elliptic_k

def solve_problem():
    """
    Solves the problem by finding the smallest n where f(n) > 10 and then
    calculating n * ||W_n||_inf.
    """
    # Step 1: Symbolically define the function and generate its Taylor series coefficients.
    # The function is g(x) = (2/pi) * K(x) * exp(x), where K(x) is the complete
    # elliptic integral of the first kind. In sympy, K(m) is elliptic_k(m).
    x = sp.Symbol('x')
    g_expr = (2 / sp.pi) * elliptic_k(x) * sp.exp(x)

    # We will need to test n up to a reasonable limit. Let's pre-compute
    # coefficients for n up to 30.
    max_degree = 30
    taylor_coeffs_sympy = []
    # The series needs to be computed up to order n to get the nth coefficient.
    series_expr = g_expr.series(x, 0, max_degree + 1)
    for i in range(max_degree + 1):
        coeff = series_expr.coeff(x, i)
        taylor_coeffs_sympy.append(coeff)

    # Convert the symbolic coefficients (fractions) to floating-point numbers.
    taylor_coeffs = [float(c) for c in taylor_coeffs_sympy]

    # Step 2: Iterate through n to find the smallest n where f(n) > 10.
    for n in range(1, max_degree + 1):
        # The Taylor polynomial P_n(x) has degree n, so it includes coefficients from a_0 to a_n.
        poly_coeffs_for_n = taylor_coeffs[0:n + 1]

        # numpy.roots requires coefficients from the highest power of x to the lowest.
        # Our list is [a_0, a_1, ..., a_n], so we need to reverse it.
        reversed_coeffs = list(reversed(poly_coeffs_for_n))

        # The eigenvalues of the companion matrix are the roots of the polynomial.
        eigenvalues = np.roots(reversed_coeffs)

        # Step 3: Calculate f(n), the sum of the absolute cubes of the eigenvalues.
        f_n = np.sum(np.abs(eigenvalues)**3)

        # Step 4: Check if the condition f(n) > 10 is met.
        if f_n > 10:
            target_n = n
            
            # Step 5: Calculate the infinity norm of W_n.
            # Assuming distinct eigenvalues, W_n is a diagonal matrix, so its
            # infinity norm is the maximum absolute value of the eigenvalues (spectral radius).
            inf_norm_Wn = np.max(np.abs(eigenvalues))
            
            # Step 6: Calculate the final result.
            result = target_n * inf_norm_Wn
            
            # Print the final equation as requested.
            print(f"The smallest n such that f(n) > 10 is n = {target_n}.")
            print(f"For this n, the sum of absolute cubes of eigenvalues f(n) is: {f_n}")
            print(f"The infinity norm of W_n is ||W_n||_inf = {inf_norm_Wn}")
            print(f"The final result is n * ||W_n||_inf:")
            print(f"{target_n} * {inf_norm_Wn} = {result}")
            
            # Output the final answer in the required format.
            print(f"\n<<<{result}>>>")
            return

# Execute the solver function.
solve_problem()