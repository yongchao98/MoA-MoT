import math

def solve():
    """
    This function determines the constant factor in the upper bound for ||B Q_{0,M}||_infty.

    The problem asks for the upper-bound for ||B Q_{0, M}||_infty expressed as a factor of sqrt(N).
    This suggests a bound of the form: k * sqrt(N) * Expression.
    The task is to find the constant k.

    1. The factor sqrt(N) typically arises from the matrix norm inequality:
       ||A||_infty <= sqrt(N) * ||A||_2 for a matrix A of size N x N.

    2. This implies that the proof likely first establishes a bound for the spectral norm ||B Q_{0,M}||_2
       and then converts it to the infinity norm.

    3. Let's assume the bound on the spectral norm is of the form:
       ||B Q_{0,M}||_2 <= K_prime * product_term,
       where product_term is related to the product given in the text, i.e., product_{t=0 to M}(1 - c*delta_t).

    4. Applying the norm conversion, we get:
       ||B Q_{0,M}||_infty <= sqrt(N) * ||B Q_{0,M}||_2 <= sqrt(N) * K_prime * product_term.

    5. The constant factor of sqrt(N) is therefore K_prime.

    6. In such theoretical analyses, the leading constant K_prime is often 1. This is because the recursive
       arguments often start with a normalized quantity, like the projection matrix B, which has a spectral
       norm ||B||_2 = 1. Without a full proof providing another constant, the most reasonable assumption is K_prime = 1.

    7. Therefore, the upper bound is 1 * sqrt(N) * product_term.
    """

    # The constant factor k in the expression k * sqrt(N) * ...
    k = 1

    # The final equation for the bound is ||B Q_{0, M}||_infty <= k * sqrt(N) * product.
    # The problem asks to output the factor of sqrt(N).
    print(k)

solve()