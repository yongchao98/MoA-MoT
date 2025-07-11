def compute_bounded_cohomology_dim():
    """
    Computes the dimension of the degree 4 bounded cohomology group
    of T x T, where T is Thompson's group.
    """
    # Let d_k = dim(H_b^k(T, R)), the dimension of the k-th bounded
    # cohomology group of Thompson's group T with trivial real coefficients.
    # The known dimensions are:
    # d_0 = 1
    # d_1 = 0
    # d_2 = 1
    # d_k = 0 for k >= 3
    # We use a dictionary to store the non-zero dimensions.
    # d.get(k, 0) will return the dimension, or 0 if not specified.
    d = {
        0: 1,
        1: 0,
        2: 1,
    }

    # We want to compute the dimension for degree n=4.
    n = 4

    print("To compute the dimension of the degree 4 bounded cohomology group of T x T,")
    print("we use the Künneth formula for bounded cohomology.")
    print("Let d_k = dim H_b^k(T, R). The formula is:")
    print(f"Dimension = d_0*d_4 + d_1*d_3 + d_2*d_2 + d_3*d_1 + d_4*d_0\n")
    
    total_dim = 0
    equation_parts = []

    # The sum is over p from 0 to n: sum_{p=0 to n} d_p * d_{n-p}
    for p in range(n + 1):
        q = n - p
        d_p = d.get(p, 0)
        d_q = d.get(q, 0)
        term = d_p * d_q
        total_dim += term
        equation_parts.append(f"{d_p}*{d_q}")

    final_equation = " + ".join(equation_parts)
    
    print("Substituting the known values:")
    print(f"d_0={d.get(0,0)}, d_1={d.get(1,0)}, d_2={d.get(2,0)}, d_3={d.get(3,0)}, d_4={d.get(4,0)}\n")

    print("The calculation is:")
    print(f"Dimension = {final_equation}")
    print(f"          = {total_dim}")
    
    print(f"\nThe dimension of the degree 4 bounded cohomology group of T x T is {total_dim}.")

compute_bounded_cohomology_dim()