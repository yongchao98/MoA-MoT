import itertools

def solve_letter_arrangement():
    """
    Calculates the number of valid arrangements for the letters L, N, S, W.
    """
    # Step 1: Define letter properties based on tip positions.
    # (Left Tip, Right Tip), where Top=0, Bottom=1.
    letter_tips = {
        'L': (0, 1),  # Top -> Bottom
        'N': (1, 0),  # Bottom -> Top
        'S': (0, 1),  # Top -> Bottom
        'W': (0, 0)   # Top -> Top
    }

    # Step 2: Define the connection rule.
    # A connection is valid if the right tip of the first letter and the
    # left tip of the second letter are at different vertical levels.
    def can_connect(char1, char2):
        """Checks if char1 can connect to char2."""
        right_tip1 = letter_tips[char1][1]
        left_tip2 = letter_tips[char2][0]
        return right_tip1 != left_tip2

    # Step 3: Iterate through all permutations to find valid arrangements.
    letters = ['L', 'N', 'S', 'W']
    valid_arrangements = []
    all_permutations = itertools.permutations(letters)

    for p in all_permutations:
        # An arrangement is valid if each adjacent pair can connect.
        if (can_connect(p[0], p[1]) and
            can_connect(p[1], p[2]) and
            can_connect(p[2], p[3])):
            valid_arrangements.append(" -> ".join(p))

    # Step 4: Print the results and the final equation.
    print("Found the following valid arrangements:")
    for arr in valid_arrangements:
        print(arr)
    
    count = len(valid_arrangements)
    # Create the equation string, e.g., "1 + 1 = 2"
    if count > 0:
        equation_parts = ["1"] * count
        equation_str = " + ".join(equation_parts)
        print(f"\nThe final count is {equation_str} = {count}")
    else:
        print("\nNo valid arrangements found. The total is 0.")

solve_letter_arrangement()