import itertools

def solve_letter_arrangement():
    """
    Finds and counts the number of valid arrangements of the letters L, N, S, W
    based on a specific connection rule.
    """
    # Step 1: Define the properties for each letter.
    # The tuple represents (position_of_leftmost_tip, position_of_rightmost_tip).
    letter_properties = {
        'L': ('Top', 'Bottom'),
        'N': ('Bottom', 'Top'),
        'S': ('Top', 'Bottom'),
        'W': ('Top', 'Top')
    }
    
    letters = ['L', 'N', 'S', 'W']
    valid_arrangements = []
    
    # Step 2: Generate all possible arrangements (permutations) of the letters.
    all_permutations = list(itertools.permutations(letters))
    
    # Step 3: Iterate through each permutation to check its validity.
    for p in all_permutations:
        is_arrangement_valid = True
        # An arrangement is a sequence like ('L', 'S', 'W', 'N').
        # We need to check the connection for each adjacent pair.
        # e.g., L->S, S->W, W->N
        for i in range(len(p) - 1):
            letter1 = p[i]
            letter2 = p[i+1]
            
            right_tip_1 = letter_properties[letter1][1]
            left_tip_2 = letter_properties[letter2][0]
            
            # Step 4: Apply the connection rule.
            # The connection is valid if the tip positions are different.
            if right_tip_1 == left_tip_2:
                is_arrangement_valid = False
                break  # This arrangement is invalid, move to the next one.
        
        if is_arrangement_valid:
            valid_arrangements.append(p)

    # Step 5: Print the results.
    print("Found the following valid arrangements:")
    for arr in valid_arrangements:
        # We print each letter in the equation that forms the arrangement
        print(f"    {' -> '.join(arr)}")
        
    print(f"\nThe total number of ways to arrange the letters is: {len(valid_arrangements)}")

solve_letter_arrangement()
<<<2>>>