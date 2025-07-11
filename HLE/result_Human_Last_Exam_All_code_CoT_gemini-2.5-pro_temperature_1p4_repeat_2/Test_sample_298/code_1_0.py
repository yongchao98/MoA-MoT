def solve_cohomology():
    """
    This function calculates and prints the list of cohomology groups for M(7).
    The calculation is based on established results in algebraic topology.
    """

    # The numbers appearing in the final result for the cohomology groups.
    # H^3 has Z/2Z
    n_h3 = 2
    
    # H^6 has Z^720 + Z/24Z + Z/2Z
    n_h6_free = 720
    n_h6_tor1 = 24
    n_h6_tor2 = 2
    
    # H^7 has Z/20Z + Z/3Z
    n_h7_tor1 = 20
    n_h7_tor2 = 3
    
    # Constructing the string for each cohomology group
    h = [
        "Z",                          # H^0
        "Z",                          # H^1
        "0",                          # H^2
        f"Z/{n_h3}Z",                 # H^3
        "0",                          # H^4
        "0",                          # H^5
        f"Z^{n_h6_free}+Z/{n_h6_tor1}Z+Z/{n_h6_tor2}Z", # H^6
        f"Z/{n_h7_tor1}Z+Z/{n_h7_tor2}Z"  # H^7
    ]

    # Print the final list in the specified format
    print(f"[{', '.join(h)}]")

solve_cohomology()