import itertools

def solve_letter_arrangement():
    """
    Finds the number of ways to arrange "L", "N", "S", "W" based on a connection rule.

    The connection rule is based on the type of the leftmost and rightmost unconnected tips
    of the capital letters.
    """
    letters = ['L', 'N', 'S', 'W']

    # Define the type of the leftmost and rightmost unconnected tips for each letter.
    # Categories are 'top', 'bottom', and 'middle'.
    tips = {
        'L': {'left': 'top', 'right': 'middle'},
        'N': {'left': 'top', 'right': 'bottom'},
        'S': {'left': 'bottom', 'right': 'top'},
        'W': {'left': 'top', 'right': 'top'}
    }

    def can_connect(letter1, letter2):
        """Checks if letter1 can connect to letter2."""
        return tips[letter1]['right'] == tips[letter2]['left']

    # Generate all permutations of the four letters.
    all_permutations = list(itertools.permutations(letters))

    valid_arrangements = []
    
    # Check each permutation for validity.
    for p in all_permutations:
        # A valid arrangement p = (p1, p2, p3, p4) must satisfy three connections.
        if can_connect(p[0], p[1]) and can_connect(p[1], p[2]) and can_connect(p[2], p[3]):
            valid_arrangements.append("".join(p))

    # Print the results
    if not valid_arrangements:
        print("No valid arrangements found.")
        print("0")
    else:
        print("The valid arrangements are:")
        for arr in valid_arrangements:
            print(arr)
        
        # Create and print the final equation as requested.
        equation_parts = ['1'] * len(valid_arrangements)
        equation_str = " + ".join(equation_parts)
        print(f"{equation_str} = {len(valid_arrangements)}")


solve_letter_arrangement()