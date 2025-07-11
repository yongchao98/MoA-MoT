def dim_Hb_T(k):
    """
    Computes the dimension of the k-th bounded cohomology group of Thompson's group T
    with trivial real coefficients.
    - dim H_b^0(T, R) = 1
    - dim H_b^k(T, R) = 0 for k >= 1
    """
    if k == 0:
        return 1
    else:
        return 0

def compute_dimension_product(group1_dim_func, group2_dim_func, n):
    """
    Computes the dimension of the n-th bounded cohomology group of a product of two groups
    using the Künneth formula.
    """
    total_dim = 0
    
    print(f"To compute dim H_b^{n}(T x T, R), we use the Künneth formula:")
    print(f"dim H_b^{n}(T x T, R) = sum_{{p+q={n}}} dim(H_b^p(T, R)) * dim(H_b^q(T, R))\n")
    
    equation_parts = []
    for p in range(n + 1):
        q = n - p
        dim_p = group1_dim_func(p)
        dim_q = group2_dim_func(q)
        term_dim = dim_p * dim_q
        
        print(f"Term p={p}, q={q}:")
        print(f"  dim(H_b^{p}(T, R)) * dim(H_b^{q}(T, R)) = {dim_p} * {dim_q} = {term_dim}")
        
        total_dim += term_dim
        equation_parts.append(str(term_dim))
        
    final_equation = " + ".join(equation_parts)
    print(f"\nThe total dimension is the sum of all terms:")
    print(f"Total Dimension = {final_equation} = {total_dim}")

    return total_dim

if __name__ == '__main__':
    # We want to compute the dimension for degree n=4
    degree_n = 4
    
    # Compute the dimension
    result = compute_dimension_product(dim_Hb_T, dim_Hb_T, degree_n)
    
    print(f"\nThe dimension of the degree {degree_n} bounded cohomology group of T x T is: {result}")
