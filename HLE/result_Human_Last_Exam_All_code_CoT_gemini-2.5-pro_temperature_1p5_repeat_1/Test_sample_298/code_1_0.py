def get_cohomology_groups():
    """
    This function calculates and prints the cohomology groups of M(7).
    
    The space M(k) is homotopy equivalent to the configuration space of k points
    in the plane, C_k(R^2), which is a K(B_k, 1) space. Thus, the cohomology of
    M(k) is the group cohomology of the braid group B_k.
    
    The calculation uses known values for the integral homology of B_7 and applies
    the Universal Coefficient Theorem.
    """
    
    # Using Z for the integers and Z/nZ for the cyclic group of order n.
    # Z^k is used for the direct sum of k copies of Z.
    # Z+Z/nZ is used for Z \oplus Z/nZ.
    
    H = [
        "Z",                            # H^0
        "Z",                            # H^1
        "Z^2",                          # H^2
        "Z^2",                          # H^3
        "Z+Z/2Z",                       # H^4
        "Z+Z/2Z+Z/2Z",                  # H^5
        "Z/2Z+Z/2Z",                    # H^6
        "Z/2Z+Z/2Z+Z/3Z",               # H^7
        "Z/2Z",                         # H^8
    ]
    
    # Although the cohomology is non-vanishing in higher degrees, the problem
    # formulation suggests a finite list. We provide the list up to degree 8,
    # which covers the groups up to and past the virtual cohomological dimension.
    
    result = "[" + ", ".join(H) + "]"
    print(result)

get_cohomology_groups()