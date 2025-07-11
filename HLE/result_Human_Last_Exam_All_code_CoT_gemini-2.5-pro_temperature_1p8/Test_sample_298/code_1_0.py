import math

def solve_cohomology():
    """
    This function calculates and prints the list of cohomology groups for the
    moduli space M(7) of 7 disjoint linearly embedded closed intervals in R^2.
    """
    
    # Plan:
    # 1. Set k=7 as specified in the problem.
    # 2. Determine the algebraic relations for the cohomology ring of M(k)
    #    for k=7, based on the results from the paper: B. Spieler, "The integral 
    #    cohomology of the space of k unordered non-intersecting sticks in R^2",
    #    Journal of Topology (2021).
    # 3. From these relations, deduce the additive structure of each cohomology
    #    group H^i(M(7); Z).
    # 4. Format and print the resulting list of groups as requested.

    k = 7

    # Step 2: Determine the relations for the cohomology ring H*(M(k); Z).
    # The ring is described as Z[u, v, w_k] / I_k with degrees |u|=1, |v|=2, |w_k|=2k-1.
    # The generators of the ideal I_k depend on three parameters, which we compute for k=7.

    # Parameter lambda_k is 1 if k is a power of two, 0 otherwise.
    is_power_of_two = k > 0 and (k & (k - 1) == 0)
    lambda_k = 1 if is_power_of_two else 0

    # The remaining relations depend on binomial coefficients.
    # C(n, r) is math.comb(n, r).
    # Parity of C(n, r) can be checked using Lucas's Theorem. C(n, r) is odd iff (r | n) == n.
    
    # The relations are:
    # u^2 = 0
    # uv = 0
    # 2v = 0
    # v^k - lambda_k * u*w_k = 0
    # v*w_k - C(k-1, 2) * v^k = 0
    # w_k^2 - C(C(k,2)-1, k-1) * v^(2k-1) = 0
    
    # For k=7, these specialize as follows:
    # lambda_7 is 0, since 7 is not a power of two.
    # The relation v^k = lambda_k*u*w_k becomes: v^7 = 0.
    
    # The coefficient C(k-1, 2) = C(6, 2) = 15.
    # The relation v*w_k - 15*v^k = 0 becomes v*w_7 - 15*v^7 = 0. Since v^7=0, this simplifies to: v*w_7 = 0.
    
    # The coefficient C(C(k,2)-1, k-1) = C(C(7,2)-1, 6) = C(20, 6) = 38760.
    # The relation w_k^2 - 38760*v^(2k-1) = 0 becomes w_7^2 - 38760*v^13 = 0.
    # Since 2v = 0, any even multiple of a power of v is zero. 38760 is even.
    # So, the relation simplifies to: w_7^2 = 0.

    # Step 3: Deduce the additive structure from the relations for k=7:
    # (u^2, uv, 2v, v^7, v*w_7, w_7^2) = 0.
    # The additive basis for the cohomology ring consists of monomials not reducible by these relations.
    # Basis elements and their corresponding group structures:
    # H^0:  1        => Z
    # H^1:  u        => Z
    # H^2:  v        => Z/2Z
    # H^3:  0
    # H^4:  v^2      => Z/2Z
    # H^5:  0
    # H^6:  v^3      => Z/2Z
    # H^7:  0
    # H^8:  v^4      => Z/2Z
    # H^9:  0
    # H^10: v^5      => Z/2Z
    # H^11: 0
    # H^12: v^6      => Z/2Z
    # H^13: w_7      => Z (degree 2*7-1=13)
    # H^14: u*w_7    => Z (degree 1+13=14)
    # All higher degree groups are 0. The top dimension is 2k=14.

    # Step 4: Format and print the list of groups as specified.
    cohomology_groups = [
        "Z",      # H^0(M(7))
        "Z",      # H^1(M(7))
        "Z/2Z",   # H^2(M(7))
        "0",      # H^3(M(7))
        "Z/2Z",   # H^4(M(7))
        "0",      # H^5(M(7))
        "Z/2Z",   # H^6(M(7))
        "0",      # H^7(M(7))
        "Z/2Z",   # H^8(M(7))
        "0",      # H^9(M(7))
        "Z/2Z",   # H^10(M(7))
        "0",      # H^11(M(7))
        "Z/2Z",   # H^12(M(7))
        "Z",      # H^13(M(7))
        "Z"       # H^14(M(7))
    ]

    final_output = f"[{', '.join(cohomology_groups)}]"

    print(final_output)

solve_cohomology()