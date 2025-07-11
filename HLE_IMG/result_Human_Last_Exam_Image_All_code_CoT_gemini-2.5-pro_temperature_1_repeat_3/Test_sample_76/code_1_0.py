import itertools

def solve_cubane_substitution():
    """
    Determines the four possible pairs of carbon atoms on the cubane product
    that can be substituted with a carboxylic acid, based on the provided
    double Favorskii rearrangement mechanism.
    """
    # From the analysis of the first Favorskii rearrangement (top of the molecule),
    # the involved carbons are C2 (with Br) and C6 (with H).
    # The resulting COOH group can theoretically be on either C2 or C6.
    site1_possibilities = [2, 6]

    # From the analysis of the second (quasi) Favorskii rearrangement (bottom),
    # the involved carbons are C8 (with Br) and C4 (with the likely acidic H).
    # The resulting COOH group can theoretically be on either C4 or C8.
    site2_possibilities = [4, 8]

    # The four theoretical products result from combining the possibilities from each site.
    # We use itertools.product to get all combinations.
    all_combinations = list(itertools.product(site1_possibilities, site2_possibilities))

    # Sort the pairs for consistent ordering.
    # First sort the numbers within each tuple, then sort the list of tuples.
    sorted_pairs = sorted([tuple(sorted(pair)) for pair in all_combinations])
    
    # Format the output string as requested: (a,b), (c,d), (e,f), (g,h)
    # The problem description has a typo in the format example "(f,h)" which should be "(g,h)".
    # We will follow the pattern.
    output_string = ", ".join([f"({p[0]},{p[1]})" for p in sorted_pairs])
    
    print(output_string)

solve_cubane_substitution()