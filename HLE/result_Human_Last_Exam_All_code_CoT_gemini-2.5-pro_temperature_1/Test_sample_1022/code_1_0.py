def solve_cohomology_dimension():
    """
    Computes the dimension of the degree 4 bounded cohomology group
    of T x T with trivial real coefficients.
    """
    target_degree = 4

    # The dimension of the bounded cohomology of Thompson's group T with trivial
    # real coefficients, dim(H_b^n(T; R)), is a known result in geometric group theory.
    # We store these dimensions for n=0, 1, 2, 3. For n>=4, the dimension is 0.
    dim_Hb_T = {0: 1, 1: 0, 2: 1, 3: 1}

    print(f"To compute the dimension of the degree {target_degree} bounded cohomology group of T x T,")
    print("we use the Künneth formula for bounded cohomology:")
    print(f"dim(H_b^{target_degree}(T x T)) = sum_{{p+q={target_degree}}} dim(H_b^p(T)) * dim(H_b^q(T))")
    print("-" * 30)

    total_dim = 0
    equation_parts = []
    numeric_parts = []
    term_values = []

    # Iterate through all pairs (p, q) such that p + q = target_degree
    for p in range(target_degree + 1):
        q = target_degree - p
        
        # Get dimensions, defaulting to 0 for any degree not in the dictionary (i.e., >= 4)
        dim_p = dim_Hb_T.get(p, 0)
        dim_q = dim_Hb_T.get(q, 0)
        
        term = dim_p * dim_q
        total_dim += term
        
        # Store parts of the equation for a detailed printout
        numeric_parts.append(f"({dim_p} * {dim_q})")
        term_values.append(str(term))

    print("The calculation is:")
    print(" + ".join(numeric_parts))
    print("= " + " + ".join(term_values))
    print(f"= {total_dim}")
    print("-" * 30)
    print(f"The final dimension is: {total_dim}")

solve_cohomology_dimension()
<<<1>>>