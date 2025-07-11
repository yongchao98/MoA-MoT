import itertools

def solve_letter_arrangement():
    """
    Finds the number of ways to arrange the letters L, N, S, W
    based on a tip-to-tip connection rule.
    """
    # Step 1: Define letter properties (left_tip, right_tip)
    letter_properties = {
        'L': ('top', 'bottom'),
        'N': ('bottom', 'top'),
        'S': ('top', 'bottom'),
        'W': ('top', 'top')
    }
    letters = list(letter_properties.keys())

    valid_arrangements = []
    
    # Step 2: Generate all permutations of the letters
    all_permutations = itertools.permutations(letters)

    # Step 3: Iterate through permutations and check for validity
    for p in all_permutations:
        is_valid = True
        # Check the 3 connections in the 4-letter arrangement
        for i in range(len(p) - 1):
            prev_letter = p[i]
            next_letter = p[i+1]

            # Rule: The right tip of the previous letter must match the left tip of the next one.
            if letter_properties[prev_letter][1] != letter_properties[next_letter][0]:
                is_valid = False
                break
        
        if is_valid:
            valid_arrangements.append(p)

    # Step 4: Output the results
    print("Found the following valid arrangements:")
    for arr in valid_arrangements:
        print(" -> ".join(arr))
    
    count = len(valid_arrangements)
    # Fulfilling the request to "output each number in the final equation"
    if count > 0:
        equation_str = " + ".join(['1'] * count)
        print(f"\nThe final calculation is: {equation_str} = {count}")
    
    print(f"\nThe total number of ways is: {count}")

solve_letter_arrangement()