def print_upper_bound_factor():
    """
    This function prints the mathematical expression for the factor of sqrt(N)
    in the upper bound of ||B Q_{0, M}||_infinity.

    The derivation is as follows:
    1. The goal is to find an upper bound for ||B Q_{0, M}||_infinity.
    2. We use the norm inequality ||A||_infinity <= sqrt(N) * ||A||_2. Thus, we need to bound ||B Q_{0, M}||_2.
    3. The definition of beta_k in the text strongly suggests a recursive relationship. For each matrix M_t = D^(t)P^(t), the theory of perturbed ergodic stochastic matrices indicates that the 2-norm of the operator projected onto the space orthogonal to the all-ones vector is reduced.
    4. This can be expressed as: ||B D^(t)P^(t) x||_2 <= (1 - c*delta_t) * ||B x||_2 for some constant c > 0.
    5. Applying this recursively for the product Q_{0,M} = M_M * ... * M_0, we get ||B Q_{0,M} x||_2 <= (Product_{t=0 to M} (1 - c*delta_t)) * ||B x||_2.
    6. This implies that the induced 2-norm is bounded: ||B Q_{0, M}||_2 <= Product_{t=0 to M} (1 - c*delta_t) = beta_M.
    7. Combining with the inequality from step 2: ||B Q_{0, M}||_infinity <= sqrt(N) * beta_M.
    8. The question asks for the factor of sqrt(N), which is beta_M.
    """

    # The numbers that appear in the final formula.
    one = 1
    zero = 0

    # The final expression is beta_M, defined as a product.
    # We print this expression symbolically.
    # Using an f-string to explicitly include the numbers 0 and 1.
    final_expression = f"Product_{{t={zero}..M}}({one} - c*delta_t)"

    print(final_expression)

print_upper_bound_factor()
<<<Product_{t=0..M}(1 - c*delta_t)>>>