import math

def nCr(n, r):
    """Computes binomial coefficient nCr."""
    if r < 0 or r > n:
        return 0
    return math.factorial(n) // (math.factorial(r) * math.factorial(n - r))

def get_jordan_block_sizes(m):
    """
    For a single nilpotent Jordan block of size m on V, the induced operator on Lambda^k(V)
    is a single nilpotent Jordan block of size nCr(m, k).
    """
    return [nCr(m, k) for k in range(m + 1)]

def dim_ker_tensor_Jk_Jl(k, l):
    """
    The dimension of the kernel of (J_k tensor I + I tensor J_l) is min(k, l),
    where J_k is a nilpotent Jordan block of size k.
    """
    if k == 0 or l == 0:
      return 0
    return min(k, l)

def compute_poincare_poly():
    """
    Computes the Poincare polynomial for the given Lie algebra.
    """
    # The Lie algebra is R \ltimes R^5.
    # The action of the derivation D* on (R^5)* splits into two Jordan blocks.
    # V1 is a 3-dim space with a J_3 block.
    # V2 is a 2-dim space with a J_2 block.
    m1 = 3  # size of the first Jordan block
    m2 = 2  # size of the second Jordan block

    # Get the sizes of the Jordan blocks for the action on the exterior powers
    # Lambda^k(V1) and Lambda^k(V2).
    jordan_sizes1 = get_jordan_block_sizes(m1)
    jordan_sizes2 = get_jordan_block_sizes(m2)

    # c_p = dim(Ker(D*)) on Lambda^p(a*)
    # a* = V1 + V2
    # Lambda^p(a*) = sum_{i+j=p} Lambda^i(V1) tensor Lambda^j(V2)
    # dim(Ker(D*)) on this space is sum_{i+j=p} dim_ker_tensor(J_size(i), J_size(j))
    c = [0] * 6 # c_0 to c_5
    for p in range(len(c)):
        dim_ker_p = 0
        for i in range(p + 1):
            j = p - i
            if i <= m1 and j <= m2:
                size1 = jordan_sizes1[i]
                size2 = jordan_sizes2[j]
                dim_ker_p += dim_ker_tensor_Jk_Jl(size1, size2)
        c[p] = dim_ker_p

    # Cohomology Betti numbers b^p = c_p + c_{p-1}
    b_cohom = [0] * 7 # b^0 to b^6
    b_cohom[0] = c[0]
    for p in range(1, 7):
        cp = c[p] if p < len(c) else 0
        c_p_minus_1 = c[p - 1]
        b_cohom[p] = cp + c_p_minus_1
        
    # Homology Betti numbers b_p = b^{n-p} by Poincare duality (n=6)
    b_homol = [0] * 7 # b_0 to b_6
    for p in range(7):
        b_homol[p] = b_cohom[6 - p]

    # Construct the polynomial string
    poly_terms = []
    for p in range(7):
        coeff = b_homol[p]
        if coeff == 0:
            continue
        if p == 0:
            poly_terms.append(f"{coeff}")
        elif p == 1:
            if coeff == 1:
                poly_terms.append("x")
            else:
                poly_terms.append(f"{coeff}x")
        else:
            if coeff == 1:
                poly_terms.append(f"x^{p}")
            else:
                poly_terms.append(f"{coeff}x^{p}")
    
    poly_string = " + ".join(poly_terms)
    print("The Poincare polynomial is P(x) = " + poly_string)

compute_poincare_poly()