def get_dim_hb_T(k):
    """
    Returns the dimension of the k-th bounded cohomology group of Thompson's group T
    with trivial real coefficients, dim H_b^k(T; R).
    
    A known result is that dim H_b^k(T; R) is 1 for k=0 and 0 for k > 0.
    """
    if k == 0:
        return 1
    elif k > 0:
        return 0
    else:
        return 0 # Cohomology is 0 for negative degrees

def compute_dim_product_cohomology(n):
    """
    Computes the dimension of the n-th bounded cohomology group of T x T
    using the Künneth formula. It also prints the full equation.
    """
    total_dim = 0
    calculation_parts = []

    for p in range(n + 1):
        q = n - p
        
        # Get dimensions for H_b^p(T) and H_b^q(T)
        dim_p = get_dim_hb_T(p)
        dim_q = get_dim_hb_T(q)
        
        term_dim = dim_p * dim_q
        total_dim += term_dim
        
        # Store the numeric part for the final equation string
        calculation_parts.append(f"{dim_p}*{dim_q}")

    # Print the explanation and the final equation with numbers
    print("Using the Künneth formula for dimensions of bounded cohomology:")
    print(f"dim H_b^{n}(T x T) = sum_{{p+q={n}}} (dim H_b^p(T) * dim H_b^q(T))")
    print("\nWhich expands to the following equation with the specific dimensions substituted:")
    
    equation_str = " + ".join(calculation_parts)
    print(f"dim H_b^{n}(T x T) = {equation_str}")
    
    # Print the final result
    print(f"\nThe total sum is: {total_dim}")
    print("\nSo, the dimension of the degree 4 bounded cohomology group of T x T is:")
    print(total_dim)


if __name__ == "__main__":
    degree_n = 4
    compute_dim_product_cohomology(degree_n)
