def solve():
    """
    This function calculates the constant factor for the upper-bound of ||B Q_{0, M}||_oo.

    Based on the provided text and established results in the literature,
    the upper bound for ||B Q_{0, M}||_oo can be expressed in the form C * sqrt(N).
    The problem is to find the constant C.

    Our derivation shows that ||B Q_{0, M}||_oo <= (2 * beta_{M-1} + delta_M) * sqrt(N).
    To find a constant upper-bound factor, we find the maximum possible value of the
    expression (2*beta_{M-1} + delta_M).

    - beta_{M-1} is a product of terms less than or equal to 1, so beta_{M-1} <= 1.
    - delta_M is a norm of a matrix whose entries are bounded, such that delta_M <= 1.

    Thus, the factor is bounded by 2*1 + 1 = 3.
    """
    
    # The derived constant factor for the upper bound.
    C = 3
    
    # The final expression for the bound is C * sqrt(N). The question asks for the factor C.
    print(C)

solve()