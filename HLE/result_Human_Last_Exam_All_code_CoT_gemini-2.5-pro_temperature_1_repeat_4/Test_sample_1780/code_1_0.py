import itertools

def solve_letter_arrangement():
    """
    Finds the number of ways to arrange the letters L, N, S, W based on a connection rule.
    """
    # Step 1: Define the properties of each letter's connection points (tips).
    # The vertical position ('top' or 'bottom') of the leftmost and rightmost tip is stored.
    letter_tips = {
        'L': {'left': 'top', 'right': 'bottom'},
        'N': {'left': 'bottom', 'right': 'top'},
        'S': {'left': 'top', 'right': 'bottom'},
        'W': {'left': 'top', 'right': 'top'}
    }
    
    letters = ['L', 'N', 'S', 'W']
    valid_arrangements = []
    
    # Step 2: Generate all possible arrangements (permutations).
    all_permutations = list(itertools.permutations(letters))
    
    # Step 3: Iterate through each arrangement and check if it's valid.
    for p in all_permutations:
        is_valid_arrangement = True
        # An arrangement (p1, p2, p3, p4) is valid if all sequential connections
        # (p1->p2, p2->p3, p3->p4) are valid.
        for i in range(len(p) - 1):
            first_letter = p[i]
            second_letter = p[i+1]
            
            # Step 4: Apply the connection rule.
            # The rule is that the vertical position of the right tip of the first letter
            # must be different from the vertical position of the left tip of the second letter.
            right_tip_pos = letter_tips[first_letter]['right']
            left_tip_pos = letter_tips[second_letter]['left']
            
            if right_tip_pos == left_tip_pos:
                is_valid_arrangement = False
                break # This connection is invalid, so the whole arrangement is invalid.
        
        if is_valid_arrangement:
            valid_arrangements.append(p)

    # Step 5: Output the results.
    # To satisfy the "output each number in the final equation" instruction,
    # we first show each individual arrangement that contributes to the total count.
    
    if not valid_arrangements:
        print("No valid arrangements found.")
    
    count = 0
    for arr in valid_arrangements:
        # Each valid arrangement contributes 1 to the final sum.
        print(f"Found valid arrangement: {', '.join(arr)}. This adds 1 to the total.")
        count += 1
        
    print(f"The final sum is: {count}")


solve_letter_arrangement()