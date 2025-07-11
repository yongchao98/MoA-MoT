def print_obstruction_groups():
    """
    Prints the list of group types that determine the homotopy obstructions.
    The parameters n and k are taken from the problem description.
    """
    
    print("The homotopy-theoretic obstructions for the two paths to be homotopic are classified by a group whose structure is determined by contributions from the following list of groups:")
    
    # 1. Homology of the space X
    group1 = "H_{n-1}(X)"
    
    # 2. Homotopy groups of the special orthogonal group SO(2k)
    # We need the groups pi_q for q from 1 to n.
    group2_list = "pi_q(SO(2k)) for q = 1, 2, ..., n"
    
    print(f"\n1. {group1}")
    print(f"2. {group2_list}")

print_obstruction_groups()