import itertools

def solve_letter_arrangement():
    """
    Calculates the number of valid arrangements of the letters L, N, S, W
    based on a connection rule derived from the problem description.
    """
    # Step 1: Define the properties of each letter's unconnected tips (Top/Bottom).
    properties = {
        'L': {'left': 'T', 'right': 'B'},
        'N': {'left': 'B', 'right': 'T'},
        'S': {'left': 'T', 'right': 'B'},
        'W': {'left': 'T', 'right': 'T'}
    }
    letters = ['L', 'N', 'S', 'W']

    # Step 2: Define the connection rule based on the example (S->W is valid, W->S is not).
    def is_valid_connection(letter1, letter2):
        """Checks if letter1 can connect to letter2."""
        right_tip1 = properties[letter1]['right']
        left_tip2 = properties[letter2]['left']
        
        # A Top-to-Top connection is invalid. All others are valid.
        if right_tip1 == 'T' and left_tip2 == 'T':
            return False
        return True

    # Step 3 & 4: Generate all permutations and check their validity.
    valid_arrangements = []
    all_permutations = itertools.permutations(letters)

    for p in all_permutations:
        # A permutation p is a tuple, e.g., ('L', 'S', 'W', 'N')
        is_fully_valid = True
        # Check the 3 connections in the 4-letter sequence
        for i in range(len(p) - 1):
            if not is_valid_connection(p[i], p[i+1]):
                is_fully_valid = False
                break
        
        if is_fully_valid:
            valid_arrangements.append(p)
    
    # Calculate the breakdown and the final equation.
    counts = {letter: 0 for letter in letters}
    for arr in valid_arrangements:
        counts[arr[0]] += 1
    
    total_count = len(valid_arrangements)
    
    # Print the "equation" by showing the contribution of each starting letter.
    count_l = counts['L']
    count_n = counts['N']
    count_s = counts['S']
    count_w = counts['W']

    print(f"Number of arrangements starting with L: {count_l}")
    print(f"Number of arrangements starting with N: {count_n}")
    print(f"Number of arrangements starting with S: {count_s}")
    print(f"Number of arrangements starting with W: {count_w}")
    print("\nThe final equation for the total number of arrangements is:")
    print(f"{count_l} + {count_n} + {count_s} + {count_w} = {total_count}")
    
solve_letter_arrangement()