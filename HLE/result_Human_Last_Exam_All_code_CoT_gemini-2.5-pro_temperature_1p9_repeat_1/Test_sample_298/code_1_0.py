import sys

def get_cohomology_of_M(k):
    """
    Computes the list of cohomology groups for M(k) based on established
    topological equivalences and results.

    The moduli space M(k) of k disjoint line segments (up to SE(2) action) is
    homotopy equivalent to the configuration space of k points on a circle, C_k(S^1).
    
    The cohomology of C_k(S^1) for k>=3 is known to be:
    H^0 = Z
    H^1 = Z
    H^2 = Z/kZ
    H^i = 0 for i > 2
    """
    
    if k < 3:
        # The case k=7 is specifically asked, so we don't need to handle k<3.
        # However, for completeness, these cases are different.
        print(f"The general formula used here is valid for k>=3. The case k={k} is not handled.", file=sys.stderr)
        return []
        
    # H^0(M(k)) is Z
    h0 = "Z"
    
    # H^1(M(k)) is Z
    h1 = "Z"
    
    # H^2(M(k)) is Z/kZ
    h2 = f"Z/{k}Z"
    
    # Higher cohomology groups are trivial, so the list ends here.
    return [h0, h1, h2]

def main():
    """
    Main function to solve the problem for k=7.
    """
    k = 7
    
    groups = get_cohomology_of_M(k)
    
    # Format the output string as per the user's request.
    # The list format is [H^0, H^1, ..., H^a]
    output_string = "[" + ", ".join(groups) + "]"
    
    print(f"The list of cohomology groups of M({k}) is:")
    print(output_string)

if __name__ == "__main__":
    main()
