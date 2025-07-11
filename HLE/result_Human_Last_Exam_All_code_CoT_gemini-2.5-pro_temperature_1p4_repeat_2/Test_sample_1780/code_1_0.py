import itertools

def count_arrangements():
    """
    Calculates the number of valid arrangements of the letters L, N, S, W
    based on a specific tip-connection rule.
    """
    # Step 1: Define the leftmost and rightmost tip types for each letter.
    # The tips are classified as 'top' or 'bottom'.
    tips = {
        'L': ('top', 'bottom'),  # (leftmost_tip, rightmost_tip)
        'N': ('bottom', 'top'),
        'S': ('top', 'bottom'),
        'W': ('top', 'top')
    }

    def can_connect(letter1, letter2):
        """
        Step 2: Check if letter1 can connect to letter2 based on the rule
        derived from the problem's example ("S can connect to W" but not vice versa).
        The rule is that the connecting tips must be of different types (top != bottom).
        """
        rightmost_tip1 = tips[letter1][1]
        leftmost_tip2 = tips[letter2][0]
        return rightmost_tip1 != leftmost_tip2

    letters = ['L', 'N', 'S', 'W']
    valid_arrangement_count = 0

    # Step 3: Generate all 4! = 24 permutations of the letters.
    all_permutations = itertools.permutations(letters)

    # Step 4: Iterate through each permutation and check if it's a valid chain.
    for perm in all_permutations:
        # A valid arrangement requires all adjacent letters to connect correctly.
        # perm = (c1, c2, c3, c4)
        if can_connect(perm[0], perm[1]) and \
           can_connect(perm[1], perm[2]) and \
           can_connect(perm[2], perm[3]):
            valid_arrangement_count += 1
            
    # The final result is the total count of valid arrangements.
    print(valid_arrangement_count)

count_arrangements()
<<<2>>>