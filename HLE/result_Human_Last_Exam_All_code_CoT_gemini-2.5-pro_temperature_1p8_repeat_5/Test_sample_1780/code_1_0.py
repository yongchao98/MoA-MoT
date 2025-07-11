import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways to arrange the letters "L", "N", "S", "W"
    based on a specific connection rule.
    """
    # Step 1: Define the tip properties for each letter.
    # (leftmost_tip, rightmost_tip) where tips can be 'Top' or 'Bottom'.
    tips = {
        'L': ('Top', 'Bottom'),
        'N': ('Bottom', 'Top'),
        'S': ('Top', 'Bottom'),
        'W': ('Top', 'Top')
    }

    letters = ['L', 'N', 'S', 'W']
    valid_arrangements = []

    # Step 2: Generate all possible permutations of the four letters.
    all_permutations = itertools.permutations(letters)

    # Step 3: Iterate through each permutation and check if it's valid.
    for p in all_permutations:
        # A permutation (p0, p1, p2, p4) is valid if each adjacent pair connects.
        # The connection rule derived from the example is that the tip types must be different.
        
        # Check connection from p[0] to p[1]
        is_valid1 = tips[p[0]][1] != tips[p[1]][0]
        
        # Check connection from p[1] to p[2]
        is_valid2 = tips[p[1]][1] != tips[p[2]][0]
        
        # Check connection from p[2] to p[3]
        is_valid3 = tips[p[2]][1] != tips[p[3]][0]

        if is_valid1 and is_valid2 and is_valid3:
            valid_arrangements.append("".join(p))
    
    # Step 4: Print the final count of valid arrangements.
    # The valid arrangements found are: ['LSWN', 'SLWN']
    print(len(valid_arrangements))

solve_letter_arrangement()