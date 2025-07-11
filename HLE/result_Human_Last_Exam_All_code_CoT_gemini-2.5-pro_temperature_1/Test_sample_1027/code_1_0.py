def dim_H_BS12(k):
    """
    Returns the dimension of the k-th homology group of BS(1,2) with real coefficients.
    """
    if k == 0:
        return 1
    if k == 1:
        return 1
    return 0

def dim_H_BS21(k):
    """
    Returns the dimension of the k-th homology group of BS(2,1) with real coefficients.
    """
    if k == 0:
        return 1
    return 0

def compute_dim_homology_product(n):
    """
    Computes the dimension of the n-th homology group of G = BS(1,2) x BS(2,1)
    using the Kunneth formula.
    """
    total_dim = 0
    print(f"Computing dim H_{n}(G, R) using the Kunneth formula:")
    print(f"dim H_{n}(G, R) = sum_{{p+q={n}}} dim H_p(BS(1,2), R) * dim H_q(BS(2,1), R)")
    print("-" * 20)

    for p in range(n + 1):
        q = n - p
        dim_p = dim_H_BS12(p)
        dim_q = dim_H_BS21(q)
        term_dim = dim_p * dim_q
        
        # We only print the terms that could potentially be non-zero
        if dim_p > 0 or dim_q > 0:
            print(f"p={p}, q={q}: dim H_{p}(BS(1,2)) * dim H_{q}(BS(2,1)) = {dim_p} * {dim_q} = {term_dim}")

        total_dim += term_dim
    
    print("-" * 20)
    print(f"The total sum is: {total_dim}")
    return total_dim

# Compute the dimension for n=31
n = 31
result = compute_dim_homology_product(n)
print("\nFinal Result:")
print(f"The dimension of the homology of G with trivial real coefficients in degree {n} is {result}.")
