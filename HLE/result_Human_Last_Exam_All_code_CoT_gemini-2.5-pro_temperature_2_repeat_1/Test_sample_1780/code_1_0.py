import itertools

def solve_letter_arrangement():
    """
    Finds and counts the number of valid arrangements for the letters L, N, S, W
    based on a specific connection rule.
    """
    letters = ['L', 'N', 'S', 'W']

    # Step 1 & 2: Define letter tip properties and the connection rule.
    # (Left Tip Position, Right Tip Position) using 'T' for TOP and 'B' for BOTTOM.
    tip_positions = {
        'L': ('T', 'B'),
        'N': ('B', 'T'),
        'S': ('T', 'B'),
        'W': ('T', 'T')
    }

    def can_connect(char1, char2):
        """
        Checks if char1 can connect to char2.
        The rule is that the vertical positions of the connecting tips must be different.
        """
        right_tip_char1 = tip_positions[char1][1]
        left_tip_char2 = tip_positions[char2][0]
        return right_tip_char1 != left_tip_char2

    # Step 3: Find all valid arrangements by checking all permutations.
    valid_arrangements = []
    all_permutations = itertools.permutations(letters)

    for p in all_permutations:
        # A permutation (c1, c2, c3, c4) is valid if all connections are valid.
        is_valid_sequence = all([
            can_connect(p[0], p[1]),
            can_connect(p[1], p[2]),
            can_connect(p[2], p[3])
        ])

        if is_valid_sequence:
            valid_arrangements.append("".join(p))
            
    # Step 4: Format the output as an equation showing the sum.
    count = len(valid_arrangements)

    if count > 0:
        # Creates a string like "1 + 1" for each valid arrangement.
        sum_str = ' + '.join(['1'] * count)
        # To satisfy the instruction to "output each number in the final equation",
        # we present the count as the result of a sum.
        print(f"{sum_str} = {count}")
    else:
        print(0)

solve_letter_arrangement()