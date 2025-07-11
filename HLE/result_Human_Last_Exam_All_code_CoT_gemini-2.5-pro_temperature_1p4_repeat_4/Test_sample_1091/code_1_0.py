import sympy as sp

def solve_limit_probability():
    """
    This function analytically calculates the limit of n*P(n) as n->inf
    using the Central Limit Theorem and symbolic mathematics.
    """
    # The total number of vectors is n = 6k. For the limit n->inf, we take k->inf.
    # We will use k as a symbolic variable to represent a large integer.
    k = sp.Symbol('k', positive=True, integer=True)
    n = 6 * k

    # According to the multidimensional Central Limit Theorem, the sum vector S
    # is approximately Gaussian with mean 0 and covariance matrix Sigma.
    # Sigma is the sum of the covariance matrices of individual random vectors (epsilon_i * v_i).
    # For a single term, the covariance is E[(epsilon_i*v_i)*(epsilon_i*v_i)^T] = v_i * v_i^T.
    # We sum these over all n vectors.

    print("This script uses the Central Limit Theorem to find the limit.")
    print("1. Calculate the covariance matrix of the sum S.")
    print("2. Determine the resulting 2D Gaussian PDF.")
    print("3. Integrate the PDF over the region ||S||^2 <= 2 to find P(n).")
    print("4. Compute the limit of n*P(n) as n -> infinity.\n")

    # There are 2k vectors of each of the three types.
    num_vectors_per_type = 2 * k

    # Define the vectors as sympy Matrices
    v_A = sp.Matrix([1, 0])
    v_B = sp.Matrix([sp.S(1)/2, sp.sqrt(3)/2])
    v_C = sp.Matrix([-sp.S(1)/2, sp.sqrt(3)/2])

    # The total covariance matrix Sigma is the sum of individual covariance matrices.
    # Sigma = 2k * (v_A*v_A^T + v_B*v_B^T + v_C*v_C^T)
    Sigma_sum_term = v_A * v_A.T + v_B * v_B.T + v_C * v_C.T
    
    Sigma = num_vectors_per_type * Sigma_sum_term

    # Extract variance and covariance terms
    var_sx = Sigma[0, 0]
    var_sy = Sigma[1, 1]
    cov_sx_sy = Sigma[0, 1]

    print("Step 1: Covariance Matrix Calculation")
    print(f"Var(S_x) = {var_sx}")
    print(f"Var(S_y) = {var_sy}")
    print(f"Cov(S_x, S_y) = {cov_sx_sy}")
    print(f"The covariance matrix is diagonal: Sigma = [[3*k, 0], [0, 3*k]].\n")

    # P(n) is the integral of the Gaussian PDF over the disk ||S||^2 <= 2.
    # PDF(s_x, s_y) = 1/(2*pi*det(Sigma)^0.5) * exp(-0.5 * S^T*Sigma^-1*S)
    #             = 1/(2*pi*3*k) * exp(-(s_x^2+s_y^2)/(2*3*k))
    # We integrate this in polar coordinates (r, theta) over a disk of radius sqrt(2):
    # P(n) = Integral from 0 to 2pi d(theta) * Integral from 0 to sqrt(2) [r * PDF] dr
    
    print("Step 2: Probability P(n) Calculation")
    r, t = sp.symbols('r t', positive=True)
    # The integral over theta gives 2*pi.
    # P(n) = 2*pi * integral_0^sqrt(2) [ r/(2*pi*3*k) * exp(-r^2/(6*k)) ] dr
    #      = 1/(3*k) * integral_0^sqrt(2) [ r * exp(-r^2/(6*k)) ] dr
    # Let u = r^2/(6k), du = 2r/(6k) dr = r/(3k) dr.
    # P(n) = integral_0^(2/6k) exp(-u) du = 1 - exp(-1/(3k))
    
    prob_n_expr = 1 - sp.exp(-1 / (3 * k))
    
    print(f"The probability P(n) as a function of k is: {prob_n_expr}\n")

    # Step 3: Compute the final limit.
    print("Step 3: Limit Calculation")
    # We want the limit of n * P(n) as k -> infinity.
    expression_to_limit = n * prob_n_expr
    print(f"We need to find the limit of n*P(n) = {expression_to_limit}")

    # For large k, the expression P(n) can be approximated using the Taylor series
    # expansion of exp(-x) around x=0, which is 1 - x + O(x^2).
    # Here, x = 1/(3*k), which is small for large k.
    # P(n) = 1 - (1 - 1/(3*k) + ...) ≈ 1/(3*k).
    # So, n*P(n) ≈ 6*k * (1/(3*k)).
    
    numerator_coeff = 6
    denominator_coeff = 3
    final_value = numerator_coeff / denominator_coeff
    
    print(f"\nFor large k, the expression approximates to the equation ({numerator_coeff}*k) / ({denominator_coeff}*k).")
    print(f"The numbers in this final equation are {numerator_coeff}, {denominator_coeff}, and the result {final_value}.")
    
    # Calculate the limit formally with sympy for confirmation
    final_limit = sp.limit(expression_to_limit, k, sp.oo)

    print(f"\nThe exact value of the limit is: {final_limit}")

if __name__ == "__main__":
    solve_limit_probability()