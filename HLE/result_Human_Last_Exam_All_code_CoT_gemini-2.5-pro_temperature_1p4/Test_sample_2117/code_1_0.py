def find_least_upper_bound():
    """
    This function implements the logic to find the least upper bound of the given expression.
    The reasoning is based on the properties of the Gaussian Hessenberg Decomposition
    of the specific Cayley-Menger matrix.
    """

    # 1. Analyze the matrix P.
    # The matrix P in the Gaussian Hessenberg decomposition of C_n = J - I
    # is unit lower triangular. This is because the decomposition algorithm
    # does not require pivoting for this specific matrix C_n.

    # 2. Calculate the eigenvalues of P.
    # The eigenvalues of a unit lower triangular matrix are its diagonal entries.
    lambda_p_max = 1
    lambda_p_min = 1

    # 3. Calculate E_P, the average eigenvalue gap of P.
    # The size of P is m x m where m = n + 2. The formula for E_P is
    # (lambda_max - lambda_min) / (m - 1).
    # Since the numerator is zero, E_P is zero for any positive integer n.
    E_P = 0  # (1 - 1) / (n + 1) = 0 for n >= 1

    # 4. Evaluate the full product.
    # The product is E_P * E_H * S_P * S_H.
    # The other terms (E_H, S_P, S_H) are finite.
    # Thus, the entire product is determined by E_P.
    # We can represent the other terms symbolically for the final equation.
    E_H_symbolic = "E_H(n)"
    S_P_symbolic = "S_P(n)"
    S_H_symbolic = "S_H(n)"
    product_value = 0

    # 5. Determine the least upper bound.
    # Since the product is 0 for all n, the set of values is {0}.
    # The least upper bound of {0} is 0.
    least_upper_bound = 0

    # Output the details of the "final equation" and the result.
    print(f"The equation for the product is: E_P * E_H * S_P * S_H")
    print(f"The term E_P evaluates to: ({lambda_p_max} - {lambda_p_min}) / (n+1) = {E_P}")
    print(f"The full product is: {E_P} * {E_H_symbolic} * {S_P_symbolic} * {S_H_symbolic} = {product_value}")
    print(f"\nSince the product is {product_value} for all positive integers n, the set of possible values is {{{product_value}}}.")
    print(f"The least upper bound of this set is {least_upper_bound}.")


if __name__ == '__main__':
    find_least_upper_bound()
