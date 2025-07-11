def get_dim_Hb_T(k):
    """
    Returns the dimension of the k-th bounded cohomology group of Thompson's group T
    with trivial real coefficients.
    """
    if k == 0:
        return 1
    if k == 2:
        return 1
    # For k=1 and k>=3, the dimension is 0.
    return 0

def compute_dim_product(group_name, n):
    """
    Computes the dimension of the n-th bounded cohomology group of G x G,
    where G is the group specified by group_name.
    """
    total_dim = 0
    
    # Build string representations for the calculation steps
    sum_terms_str_list = []
    sum_values_str_list = []
    product_results_list = []

    for i in range(n + 1):
        j = n - i
        dim_i = get_dim_Hb_T(i)
        dim_j = get_dim_Hb_T(j)
        term_product = dim_i * dim_j
        total_dim += term_product
        
        sum_terms_str_list.append(f"dim H_b^{i}({group_name})*dim H_b^{j}({group_name})")
        sum_values_str_list.append(f"{dim_i}*{dim_j}")
        product_results_list.append(str(term_product))

    # Print the detailed calculation
    print(f"To compute the dimension of the degree {n} bounded cohomology group of {group_name} x {group_name}, we use the Künneth formula:")
    print(f"dim H_b^{n}({group_name} x {group_name}) = " + " + ".join(sum_terms_str_list))
    print(f"                   = " + " + ".join(sum_values_str_list))
    print(f"                   = " + " + ".join(product_results_list))
    print(f"                   = {total_dim}")


# Set the group name and the degree for the computation
GROUP_NAME = "T"
DEGREE = 4

compute_dim_product(GROUP_NAME, DEGREE)