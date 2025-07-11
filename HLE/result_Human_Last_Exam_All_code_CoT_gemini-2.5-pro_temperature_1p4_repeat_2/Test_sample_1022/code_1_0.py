def get_dim_Hb_T(k):
    """
    Returns the dimension of the k-th bounded cohomology group of Thompson's group T
    with trivial real coefficients, based on established mathematical results.
    """
    if k == 0:
        return 1
    elif k == 2:
        return 1
    else:
        return 0

def compute_dim_product_cohomology(n):
    """
    Computes the dimension of the n-th bounded cohomology group of T x T
    using the Künneth formula.
    """
    total_dim = 0
    print(f"To compute the dimension of H_b^{n}(T x T; R), we use the Künneth formula:")
    print(f"dim(H_b^{n}(T x T)) = Sum_{{p+q={n}}} [ dim(H_b^p(T)) * dim(H_b^q(T)) ]")
    print("-" * 60)

    # We need to sum over all integer pairs (p, q) such that p+q=n
    sum_string_parts = []
    for p in range(n + 1):
        q = n - p
        dim_p = get_dim_Hb_T(p)
        dim_q = get_dim_Hb_T(q)
        term_value = dim_p * dim_q
        total_dim += term_value
        
        # Print the breakdown of each term in the summation
        print(f"Term for (p={p}, q={q}): dim(H_b^{p}(T)) * dim(H_b^{q}(T)) = {dim_p} * {dim_q} = {term_value}")
        sum_string_parts.append(str(term_value))

    # Print the final equation with all the numbers
    print("-" * 60)
    final_equation = f"dim(H_b^{n}(T x T)) = " + " + ".join(sum_string_parts)
    print(final_equation)
    print(f"Total dimension = {total_dim}")


if __name__ == "__main__":
    # The degree of the cohomology group we are interested in is 4.
    degree = 4
    compute_dim_product_cohomology(degree)