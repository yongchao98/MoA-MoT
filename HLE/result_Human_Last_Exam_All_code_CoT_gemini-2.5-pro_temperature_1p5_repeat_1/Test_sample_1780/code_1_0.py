import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways to arrange the letters L, N, S, and W
    based on a specific tip-to-tip connection rule.
    """
    
    # Step 1: Define letter properties based on visual analysis.
    # We use 0 for a 'Top' tip and 1 for a 'Bottom' tip.
    # Each letter is defined by a tuple: (left_tip_type, right_tip_type).
    letter_properties = {
        'L': (1, 1),  # (Bottom, Bottom)
        'N': (0, 1),  # (Top, Bottom)
        'S': (0, 1),  # (Top, Bottom)
        'W': (0, 0)   # (Top, Top)
    }
    
    letters = sorted(letter_properties.keys())
    
    # Step 2: Generate all permutations and check for validity.
    all_permutations = itertools.permutations(letters)
    
    valid_arrangements_by_start = {letter: 0 for letter in letters}

    for arrangement in all_permutations:
        is_valid = True
        # Check each connection in the sequence (e.g., L->N, N->S, S->W)
        for i in range(len(arrangement) - 1):
            letter1 = arrangement[i]
            letter2 = arrangement[i+1]
            
            right_tip_of_first = letter_properties[letter1][1]
            left_tip_of_second = letter_properties[letter2][0]
            
            # Step 3: Apply the connection rule.
            # The rule, deduced from the S->W example, is that the
            # connecting tip types must be different.
            if right_tip_of_first == left_tip_of_second:
                is_valid = False
                break
        
        if is_valid:
            start_letter = arrangement[0]
            valid_arrangements_by_start[start_letter] += 1
            
    # Step 4: Output the result as a descriptive equation.
    counts = [valid_arrangements_by_start[l] for l in letters]
    total = sum(counts)
    
    # Creates a string like "2 + 2 + 2 + 2" from the counts list
    equation_str = " + ".join(map(str, counts))
    
    print(f"{equation_str} = {total}")

solve_letter_arrangement()