def solve_cubane_substitution():
    """
    This function determines the four possible pairs of carbon atoms on the cubane product
    that can be substituted with a carboxylic acid, based on the mechanism of the
    double Favorskii rearrangement.
    """

    # Step 1 & 2: Identify the alpha-carbons for each ketone group.
    # The top ketone is adjacent to carbons 6 and 7.
    top_ketone_alpha_carbons = [6, 7]
    # The bottom ketone is adjacent to carbons 3 and 8.
    bottom_ketone_alpha_carbons = [3, 8]

    # Step 3: Generate all four theoretical combinations.
    # The Favorskii rearrangement attaches the new COOH group to one of the
    # original alpha-carbons. We pair each possibility from the first
    # rearrangement with each possibility from the second.
    possible_pairs = []
    for c1 in top_ketone_alpha_carbons:
        for c2 in bottom_ketone_alpha_carbons:
            # Sort the numbers within the pair for consistent ordering.
            pair = tuple(sorted((c1, c2)))
            possible_pairs.append(pair)

    # Sort the final list of pairs for a standardized output.
    # The set conversion removes duplicates, although none are expected here.
    # The final sort makes the order deterministic.
    sorted_pairs = sorted(list(set(possible_pairs)))

    # Step 4: Format the output string as requested.
    # The format is (a,b), (c,d), (e,f), (g,h)
    output_string = ", ".join([f"({p[0]},{p[1]})" for p in sorted_pairs])
    
    print(output_string)

solve_cubane_substitution()