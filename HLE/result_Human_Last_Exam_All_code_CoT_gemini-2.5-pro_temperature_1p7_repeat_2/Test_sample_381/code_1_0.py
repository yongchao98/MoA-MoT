import sympy as sp

def get_upper_bound():
    """
    This function calculates and prints the symbolic upper-bound for ||B Q_{0,M}||_infinity.
    The derivation is as follows:
    1.  The quantity to bound is ||B Q_{0,M}||_infinity.
    2.  We use the decomposition Q_{0,M} = P^{(M:0)} + E_M, where P^{(M:0)} is the product of
        stochastic matrices and E_M is the error term.
    3.  By triangle inequality, ||B Q_{0,M}||_inf <= ||B P^{(M:0)}||_inf + ||B E_M||_inf.
    4.  The term ||B P^{(M:0)}||_inf is bounded using the 2-norm and norm relations:
        ||B P^{(M:0)}||_inf <= sqrt(N) * ||B P^{(M:0)}||_2 <= sqrt(N) * lambda^(M+1),
        where lambda < 1 is the contraction rate.
    5.  The error term B E_M can be expressed as a sum, and its norm is bounded by:
        ||B E_M||_inf <= sqrt(N) * sum_{k=0 to M} lambda^(M-k) * delta_k.
    6.  Combining these, we get:
        ||B Q_{0,M}||_inf <= sqrt(N) * (lambda^(M+1) + sum_{k=0 to M} lambda^(M-k) * delta_k).
    7.  To find a uniform bound for all M, we bound the terms lambda^(M+1) <= 1 and
        lambda^(M-k) <= 1.
    8.  This yields the final upper-bound: sqrt(N) * (1 + sum_{t=0 to infinity} delta_t).
    9.  The condition beta > 0 ensures that the sum sum_{t=0 to infinity} delta_t is finite.
    """
    
    # Define the symbolic variables
    N = sp.Symbol('N', positive=True, integer=True)
    t = sp.Symbol('t', integer=True, nonnegative=True)
    # delta_t is represented as a function of t
    delta = sp.Function('delta')(t)

    # Define the infinite sum part of the bound
    # This sum is finite under the condition beta > 0
    infinite_sum_delta = sp.Sum(delta, (t, 0, sp.oo))

    # Construct the full upper-bound expression
    upper_bound = sp.sqrt(N) * (1 + infinite_sum_delta)
    
    # Print the equation representing the bound
    print("The upper-bound for ||B Q_{0,M}||_infinity is:")
    # Using sp.pprint for a more readable mathematical format
    sp.pprint(upper_bound, use_unicode=True)
    print("\nWhich in linear form is:")
    print(upper_bound)


if __name__ == '__main__':
    get_upper_bound()
