import math

def compute_dimension():
    """
    Computes the dimension of the degree 4 bounded cohomology group of T x T.
    """
    
    # We define a function to retrieve the known dimensions of H_b^k(T; R).
    # 'math.inf' is used to represent an infinite dimension.
    def get_dim_Hb_T(k):
        if k == 0: return 1
        if k == 1: return 0
        if k == 2: return math.inf
        if k == 3: return 1
        return 0  # For k >= 4, the dimension is 0.

    # The degree n for the cohomology group H_b^n(T x T).
    n = 4

    # Using the Künneth formula, we calculate the dimension.
    # The dimension is the sum of products of dimensions for p+q=n.
    
    # Build a string representation of the calculation for display.
    calc_terms = []
    result_terms = []
    total_dimension = 0

    print("The dimension is calculated using the Künneth formula:")
    print(f"dim H_b^{n}(T x T) = sum over p+q={n} of [dim H_b^p(T) * dim H_b^q(T)]\n")
    print("The full calculation is:")

    for p in range(n + 1):
        q = n - p
        dim_p = get_dim_Hb_T(p)
        dim_q = get_dim_Hb_T(q)
        
        # Format dimensions for printing, using 'inf' for infinity.
        dim_p_str = 'inf' if dim_p == math.inf else str(dim_p)
        dim_q_str = 'inf' if dim_q == math.inf else str(dim_q)
        calc_terms.append(f"({dim_p_str} * {dim_q_str})")
        
        # Dimension of a tensor product V tensor W is dim(V)*dim(W).
        # If either V or W is the zero vector space (dim 0), the product is zero.
        if dim_p == 0 or dim_q == 0:
            term_dim = 0
        else:
            term_dim = dim_p * dim_q
            
        total_dimension += term_dim
        
        # Format the result of each term for printing.
        term_res_str = 'inf' if term_dim == math.inf else str(term_dim)
        result_terms.append(term_res_str)
        
    # Print the step-by-step equation
    print(f"dim = {' + '.join(calc_terms)}")
    print(f"dim = {' + '.join(result_terms)}")
    
    # Determine the final result string
    final_dim_str = 'inf' if total_dimension == math.inf else str(total_dimension)
    
    print(f"\nThe resulting dimension is: {final_dim_str}")

compute_dimension()