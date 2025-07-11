import itertools

def solve_letter_arrangement():
    """
    Calculates the number of valid arrangements for the letters L, N, S, W
    based on a specific connection rule.
    """

    # Step 1: Define letter properties based on vertical tip positions.
    # 'start' is the leftmost unconnected tip, 'end' is the rightmost.
    letter_properties = {
        'L': {'start': 'top',    'end': 'middle'},
        'N': {'start': 'bottom', 'end': 'top'},
        'S': {'start': 'top',    'end': 'bottom'},
        'W': {'start': 'top',    'end': 'top'},
    }

    # The connection rule is that the end tip of the preceding letter must have a
    # different vertical position from the start tip of the succeeding letter.

    # Step 2: Iterate through all permutations and check for validity.
    letters_to_arrange = ['L', 'N', 'S', 'W']
    valid_arrangements = []

    for p in itertools.permutations(letters_to_arrange):
        # An arrangement (p[0], p[1], p[2], p[3]) is valid if all
        # connections p[i] -> p[i+1] are valid.
        is_valid_chain = True
        for i in range(len(p) - 1):
            letter_A = p[i]
            letter_B = p[i+1]
            
            end_A = letter_properties[letter_A]['end']
            start_B = letter_properties[letter_B]['start']
            
            if end_A == start_B:
                is_valid_chain = False
                break  # This connection is invalid, so the arrangement is invalid
                
        if is_valid_chain:
            valid_arrangements.append("".join(p))

    # Step 3: Print the results as an equation.
    # Each '1' in the sum represents one valid arrangement found.
    # The valid arrangements themselves are: LSWN and SLWN.
    
    count = len(valid_arrangements)
    
    if count == 0:
        print("0")
    else:
        # Fulfills the instruction: "output each number in the final equation"
        sum_equation = " + ".join(['1' for _ in valid_arrangements])
        print(f"{sum_equation} = {count}")

solve_letter_arrangement()