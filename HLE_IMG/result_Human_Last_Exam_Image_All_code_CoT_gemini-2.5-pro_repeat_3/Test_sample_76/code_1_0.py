import itertools

def solve_chemistry_problem():
    """
    This script determines the four possible substitution patterns for the
    double Favorskii rearrangement shown in the problem.
    """

    # Step 1: Identify the two pairs of carbons involved in the rearrangements.
    # Based on the structure, the first rearrangement involves carbons 2 and 3.
    # The second rearrangement involves carbons 7 and 8.
    reacting_pair_1 = {2, 3}
    reacting_pair_2 = {7, 8}

    # Step 2: In each rearrangement, the carboxylic acid can form on either carbon.
    # We find all combinations by taking one carbon from each pair.
    # The itertools.product function is perfect for this.
    possible_products = list(itertools.product(reacting_pair_1, reacting_pair_2))

    # Step 3: Format the pairs for consistent ordering and presentation.
    # We sort each pair internally, e.g., (7, 2) becomes (2, 7).
    # Then we sort the final list of unique pairs.
    formatted_pairs = sorted(list({tuple(sorted(p)) for p in possible_products}))

    # Step 4: Generate the final output string in the specified format.
    # The format is (a,b), (c,d), ...
    output_string_parts = []
    for pair in formatted_pairs:
        # Remember to output each number in the final equation!
        output_string_parts.append(f"({pair[0]},{pair[1]})")

    final_answer = ", ".join(output_string_parts)
    print(final_answer)

solve_chemistry_problem()