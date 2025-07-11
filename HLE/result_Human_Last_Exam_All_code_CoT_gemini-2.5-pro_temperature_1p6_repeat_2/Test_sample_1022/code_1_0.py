def solve_cohomology_dimension():
    """
    Computes the dimension of the degree 4 ordinary cohomology group of T x T.
    
    This function applies the Künneth formula for ordinary cohomology.
    The formula is: dim H^n(G1 x G2) = sum_{p+q=n} dim H^p(G1) * dim H^q(G2).
    
    For Thompson's group T, the dimensions of its ordinary cohomology groups with real
    coefficients are:
    dim H^0(T) = 1
    dim H^1(T) = 1
    dim H^2(T) = 1
    dim H^k(T) = 0 for k >= 3.
    """
    
    # Degree of the cohomology group we are interested in
    n = 4
    
    # Dimensions of H^p(T; R) for p = 0, 1, 2, 3, 4
    # d[p] corresponds to dim H^p(T; R)
    d = [1, 1, 1, 0, 0]
    
    total_dim = 0
    equation_terms = []
    
    # Iterate through all possible values of p from 0 to n
    for p in range(n + 1):
        q = n - p
        
        # Get dimensions for H^p(T) and H^q(T)
        dim_p = d[p] if p < len(d) else 0
        dim_q = d[q] if q < len(d) else 0
        
        term = dim_p * dim_q
        total_dim += term
        
        # Store the string representation for the final output equation
        equation_terms.append(f"{dim_p}*{dim_q}")
        
    # Join the terms to form the full equation string
    equation_str = " + ".join(equation_terms)
    
    print(f"The dimension is calculated by the sum: dim H^4(T x T) = sum_{{p+q=4}} dim H^p(T) * dim H^q(T)")
    print(f"The calculation is: {equation_str} = {total_dim}")

solve_cohomology_dimension()
