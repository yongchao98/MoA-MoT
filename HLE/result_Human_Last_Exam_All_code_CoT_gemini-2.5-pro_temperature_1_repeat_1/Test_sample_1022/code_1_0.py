def compute_bounded_cohomology_dimension():
    """
    Computes the dimension of the degree 4 bounded cohomology group of T x T.
    
    This function uses the Künneth formula for bounded cohomology and known results
    for the dimensions of the bounded cohomology groups of Thompson's group T.
    """
    
    # Known dimensions of the bounded cohomology of Thompson's group T: dim H_b^n(T, R)
    # Stored in a dictionary where the key is the degree n.
    dim_Hb_T = {
        0: 0,
        1: 0,
        2: 1,
        3: 1,
        4: 0,
    }

    # The degree of the cohomology group to compute.
    n = 4

    print(f"We compute dim H_b^{n}(T x T, R) using the Künneth formula:")
    print(f"dim H_b^{n}(T x T) = Σ_{{p+q={n}}} (dim H_b^p(T) * dim H_b^q(T))")
    print("-" * 50)

    total_dimension = 0
    numeric_terms = []

    # Iterate through all partitions of n into two non-negative integers p and q.
    for p in range(n + 1):
        q = n - p
        
        # Get dimensions for H_b^p(T) and H_b^q(T).
        # Degrees higher than 4 have dimension 0.
        dim_p = dim_Hb_T.get(p, 0)
        dim_q = dim_Hb_T.get(q, 0)
        
        term = dim_p * dim_q
        total_dimension += term
        
        # Store the numeric term for the final equation.
        numeric_terms.append(f"{dim_p}*{dim_q}")

    # Print the breakdown of the sum
    print("The terms of the sum correspond to (p, q) pairs summing to 4:")
    print(f"p=0, q=4: term = dim H_b^0(T) * dim H_b^4(T) = {dim_Hb_T.get(0,0)} * {dim_Hb_T.get(4,0)} = {dim_Hb_T.get(0,0) * dim_Hb_T.get(4,0)}")
    print(f"p=1, q=3: term = dim H_b^1(T) * dim H_b^3(T) = {dim_Hb_T.get(1,0)} * {dim_Hb_T.get(3,0)} = {dim_Hb_T.get(1,0) * dim_Hb_T.get(3,0)}")
    print(f"p=2, q=2: term = dim H_b^2(T) * dim H_b^2(T) = {dim_Hb_T.get(2,0)} * {dim_Hb_T.get(2,0)} = {dim_Hb_T.get(2,0) * dim_Hb_T.get(2,0)}")
    print(f"p=3, q=1: term = dim H_b^3(T) * dim H_b^1(T) = {dim_Hb_T.get(3,0)} * {dim_Hb_T.get(1,0)} = {dim_Hb_T.get(3,0) * dim_Hb_T.get(1,0)}")
    print(f"p=4, q=0: term = dim H_b^4(T) * dim H_b^0(T) = {dim_Hb_T.get(4,0)} * {dim_Hb_T.get(0,0)} = {dim_Hb_T.get(4,0) * dim_Hb_T.get(0,0)}")
    print("-" * 50)

    # Print the full equation with all the numbers.
    equation_str = " + ".join(numeric_terms)
    print("The full calculation is:")
    print(f"dim H_b^{n}(T x T) = {equation_str}")
    print(f"dim H_b^{n}(T x T) = {total_dimension}")
    print("-" * 50)
    
    print(f"The dimension of the degree {n} bounded cohomology group of T x T is {total_dimension}.")

compute_bounded_cohomology_dimension()