def solve_cohomology():
    """
    This function prints the list of cohomology groups of M(7).
    The calculation is based on the homotopy equivalence of M(k) with the
    configuration space B_k(RP^1) and the known structure of the homology
    of these spaces, combined with the Universal Coefficient Theorem.
    """
    
    # H^p is the p-th cohomology group of M(7)
    h0 = "Z"
    h1 = "Z^6"
    h2 = "Z^15+Z/2Z"
    h3 = "Z^20+Z/2Z"
    h4 = "Z^15+Z/2Z"
    h5 = "Z^6+Z/2Z"
    h6 = "Z+Z/2Z"

    result_list = [h0, h1, h2, h3, h4, h5, h6]
    
    print(f"[{', '.join(result_list)}]")

solve_cohomology()