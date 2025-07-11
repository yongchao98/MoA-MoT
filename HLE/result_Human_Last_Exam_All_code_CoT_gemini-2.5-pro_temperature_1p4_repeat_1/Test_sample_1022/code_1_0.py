def get_thompson_cohomology_dim(degree):
    """
    Returns the dimension of the ordinary cohomology group H^k(T, R) for Thompson's group T.
    
    It is a known set of results in group theory that:
    - H^0(T, R) is 1-dimensional.
    - H^k(T, R) is 0-dimensional for all k > 0.
    """
    if degree == 0:
        return 1
    else:
        return 0

def compute_product_cohomology_dim(n):
    """
    Computes the dimension of H^n(T x T, R) using the Künneth formula.
    """
    total_dim = 0
    
    print(f"To compute the dimension of the ordinary cohomology group H^{n}(T x T, R), we use the Künneth formula:")
    print(f"dim(H^{n}(T x T)) = Σ_{{p+q={n}}} dim(H^p(T)) * dim(H^q(T))")
    print("\nThe summation proceeds as follows:")
    
    for p in range(n + 1):
        q = n - p
        dim_p = get_thompson_cohomology_dim(p)
        dim_q = get_thompson_cohomology_dim(q)
        term_dim = dim_p * dim_q
        
        print(f"For p={p}, q={q}: The term is dim(H^{p}(T)) * dim(H^{q}(T)) = {dim_p} * {dim_q} = {term_dim}")
        total_dim += term_dim
        
    print(f"\nThe total dimension is the sum of these terms.")
    print(f"Result: dim(H^{n}(T x T, R)) = {total_dim}")
    return total_dim

if __name__ == "__main__":
    # The degree of the cohomology group to be computed.
    degree = 4
    compute_product_cohomology_dim(degree)