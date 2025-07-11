import math

def compute_bounded_cohomology_dim():
    """
    Computes the dimension of the degree 4 bounded cohomology group
    of T x T with trivial real coefficients.
    """
    degree = 4
    
    # Known dimensions of the bounded cohomology groups of Thompson's group T, H_b^k(T, R).
    # h_b^n(T) is infinite-dimensional for n >= 3, a result by Burger and Monod.
    h_b_T = {
        0: 1,
        1: 0,
        2: 1,
        3: math.inf,
        4: math.inf
    }

    total_dim = 0
    
    # We use the Künneth formula for dimensions:
    # dim H_b^n(T x T) = sum_{p+q=n} dim H_b^p(T) * dim H_b^q(T)
    
    print(f"The dimension of H_b^{degree}(T x T, R) is the sum of the following terms:")
    
    equation_parts = []
    
    for p in range(degree + 1):
        q = degree - p
        
        dim_p = h_b_T.get(p, 0)
        dim_q = h_b_T.get(q, 0)
        
        # In the context of tensor product of vector spaces, dim(V tensor W) = dim(V) * dim(W).
        # If one space has dimension 0, the tensor product is the zero space, so its dimension is 0.
        if dim_p == 0 or dim_q == 0:
            term_dim = 0
        else:
            term_dim = dim_p * dim_q
            
        total_dim += term_dim

        # For printing, represent infinity as 'inf'
        dim_p_str = 'inf' if math.isinf(dim_p) else str(dim_p)
        dim_q_str = 'inf' if math.isinf(dim_q) else str(dim_q)
        term_dim_str = 'inf' if math.isinf(term_dim) else str(term_dim)
        
        print(f"p={p}, q={q}: dim(H_b^{p}(T)) * dim(H_b^{q}(T)) = {dim_p_str} * {dim_q_str} = {term_dim_str}")
        equation_parts.append(f"{dim_p_str}*{dim_q_str}")

    total_dim_str = 'infinity' if math.isinf(total_dim) else str(total_dim)
    
    print("\nFull equation:")
    print(f"dim H_b^4(T x T) = " + " + ".join(equation_parts))
    
    print("\nFinal Result:")
    print(f"The total dimension is {total_dim_str}.")

compute_bounded_cohomology_dim()
