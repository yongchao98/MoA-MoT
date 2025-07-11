import itertools

def solve():
    """
    Finds the number of ways to arrange the capital letters "L", "N", "S", and "W"
    based on a specific connection rule.
    """

    # Step 1: Define the properties of the letter tips.
    # 'T' stands for a tip at the Top, 'B' for a tip at the Bottom.
    # The list for each letter is [leftmost_tip, rightmost_tip].
    letter_tips = {
        'L': ['T', 'B'],
        'N': ['B', 'T'],
        'S': ['T', 'B'],
        'W': ['T', 'T']
    }

    # Step 2: Define the connection function based on the deduced rule.
    # The rule is that the right tip of the first letter and the left tip
    # of the second letter must have different vertical positions.
    def can_connect(letter1, letter2):
        right_tip_of_first = letter_tips[letter1][1]
        left_tip_of_second = letter_tips[letter2][0]
        return right_tip_of_first != left_tip_of_second

    # Step 3: Iterate through all permutations of the letters and check for validity.
    letters = ['L', 'N', 'S', 'W']
    valid_arrangement_count = 0
    
    # Generate all 4! = 24 possible arrangements.
    for p in itertools.permutations(letters):
        # An arrangement is valid if each adjacent pair can connect.
        # We check the connection from p[0] to p[1], p[1] to p[2], and p[2] to p[3].
        if (can_connect(p[0], p[1]) and
            can_connect(p[1], p[2]) and
            can_connect(p[2], p[3])):
            valid_arrangement_count += 1
            # Optional: print the valid arrangement for verification.
            # print(f"Valid arrangement: {p[0]} -> {p[1]} -> {p[2]} -> {p[3]}")

    # The final answer is the total count of valid arrangements.
    print(valid_arrangement_count)

solve()