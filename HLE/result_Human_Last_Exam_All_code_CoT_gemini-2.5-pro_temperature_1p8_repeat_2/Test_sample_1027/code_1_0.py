def h_dim_Z(k):
    """
    Computes the dimension of the k-th homology group of the integers (Z)
    with trivial real coefficients (R).
    dim H_k(Z, R) is 1 for k=0, 1 and 0 for k > 1.
    """
    if k == 0 or k == 1:
        return 1
    else:
        return 0

def main():
    """
    Calculates the dimension of H_31(G, R) by computing dim H_31(Z^2, R).
    """
    n = 31
    total_dim = 0
    
    print(f"The problem reduces to computing the dimension of H_{n}(Z^2, R).")
    print("Using the Künneth formula, this dimension is given by the sum:")
    print(f"dim(H_{n}(Z^2, R)) = sum over p from 0 to {n} of [dim(H_p(Z, R)) * dim(H_{{n}-p}(Z, R))]")
    print("-" * 40)
    
    equation_parts = []
    
    # We build the sum term by term
    for p in range(n + 1):
        q = n - p
        dim_p = h_dim_Z(p)
        dim_q = h_dim_Z(q)
        term = dim_p * dim_q
        total_dim += term
        
        # Storing each term of the sum for the final output equation
        # The term is dim(H_p(Z, R)) * dim(H_q(Z, R))
        part_str = f"({dim_p} * {dim_q})"
        equation_parts.append(part_str)

    # Outputting the full equation with the value of each term
    # For n=31, at least one of p or q must be > 1, so one of dim_p or dim_q must be 0.
    # Therefore every term is 0.
    # e.g., for p=0, q=31: dim_p=1, dim_q=0. Term is 1*0=0.
    # for p=1, q=30: dim_p=1, dim_q=0. Term is 1*0=0.
    # for p=2, q=29: dim_p=0, dim_q=0. Term is 0*0=0.
    # etc.

    full_equation = " + ".join(equation_parts)
    print(f"The calculation of the sum is:\n{full_equation}")
    print(f"= {total_dim}")
    print("-" * 40)

    print(f"The dimension of the homology of G with trivial real coefficients in degree {n} is: {total_dim}")

if __name__ == "__main__":
    main()