import itertools

def solve_letter_arrangement():
    """
    Finds the number of ways the letters L, N, S, W can be arranged
    based on a specific tip-to-tip connection rule.
    """

    # 1. Define the properties of each letter's tips.
    # The tuple represents (leftmost_tip_position, rightmost_tip_position).
    # 'H' stands for High, 'L' stands for Low.
    letter_properties = {
        'L': ('H', 'L'),
        'N': ('L', 'H'),
        'S': ('L', 'H'),
        'W': ('H', 'H'),
    }

    letters = list(letter_properties.keys())
    
    # 2. Generate all possible arrangements (permutations) of the letters.
    all_permutations = itertools.permutations(letters)

    valid_arrangements = []
    
    # 3. Check each permutation against the connection rule.
    for p in all_permutations:
        # A permutation is a tuple, e.g., ('N', 'L', 'S', 'W')
        is_valid = True
        for i in range(len(p) - 1):
            current_letter = p[i]
            next_letter = p[i+1]

            # The rule: The rightmost tip of the current letter must match
            # the leftmost tip of the next letter.
            right_tip_current = letter_properties[current_letter][1]
            left_tip_next = letter_properties[next_letter][0]

            if right_tip_current != left_tip_next:
                is_valid = False
                break  # This arrangement is invalid, move to the next one.
        
        if is_valid:
            valid_arrangements.append(p)

    # 4. Output the calculation as an equation.
    # Each valid arrangement contributes 1 to the total count.
    if not valid_arrangements:
        print("0")
    else:
        # Create a list of '1's, one for each valid arrangement found.
        sum_components = ['1'] * len(valid_arrangements)
        # Format and print the equation and the final result.
        equation = " + ".join(sum_components)
        total = len(valid_arrangements)
        print(f"The total number of valid arrangements is {total}.")
        print(f"Calculation: {equation} = {total}")

solve_letter_arrangement()