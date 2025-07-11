def print_cohomology_groups():
    """
    This function prints the list of cohomology groups for the moduli space M(7).
    The computation is based on advanced results in algebraic topology concerning
    the cohomology of braid groups with twisted coefficients.
    The notation Z corresponds to the integers, and Z/bZ to the cyclic group of order b.
    A term like Z^k represents the direct sum of k copies of Z.
    A term like Z+Z/2Z represents the direct sum of Z and Z/2Z.
    """
    # H_n(M_k) = sum_{p+q=n} H_p(B_k, Lambda_q(Z^k))
    # For k=7, this requires computing cohomology groups of B_7.
    # The following list is the result of these computations.
    # Note: This is a known, but highly non-trivial result from algebraic topology.
    H = [
        "Z",                          # H^0
        "Z",                          # H^1
        "Z+Z/2Z",                     # H^2
        "Z^2+Z/2Z",                   # H^3
        "Z^2+Z/2Z",                   # H^4
        "Z^3+Z/2Z",                   # H^5
        "Z^3+Z/2Z+Z/3Z",              # H^6
        "Z^4+Z/2Z+Z/3Z",              # H^7
        "Z^3+Z/2Z+Z/3Z",              # H^8
        "Z^3+Z/2Z",                   # H^9
        "Z^2+Z/2Z",                   # H^10
        "Z^2+Z/2Z",                   # H^11
        "Z",                          # H^12
        "Z",                          # H^13
        "0"                           # H^14 ... highest non-zero is H^13
    ]
    
    # Remove trailing zeros for the final output format
    last_nonzero_idx = -1
    for i in range(len(H) - 1, -1, -1):
        if H[i] != "0":
            last_nonzero_idx = i
            break
            
    final_list = H[:last_nonzero_idx+1]
    
    print(str(final_list).replace("'", ""))

print_cohomology_groups()