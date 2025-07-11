def solve_rope_puzzle():
    """
    This script solves the magical rope puzzle by systematically checking all valid possibilities.
    It works backward from the final cut to establish conditions, then iterates through
    the possibilities of the first cut to find the correct sequence of events.
    """
    initial_length = 60
    target_piece = 15
    shrink_factor = 0.75
    solution_found = False

    # Iterate through possible lengths for the longer piece from the first cut (L1_long).
    # It must be > 30 cm to be the longer piece.
    # It must also be a multiple of 4 for its shrunken length to be an integer.
    for l1_long in range(32, initial_length, 4):
        
        l1_short = initial_length - l1_long
        l1_long_shrunk = int(l1_long * shrink_factor)
        
        # The piece used for the second cut (L2) must be the shrunken longer piece.
        # If we check the shorter piece (l1_short), it's always too small to satisfy the conditions.
        # For example, if l1_long=52, l1_short=8. If l1_long=48, l1_short=12. Not enough to cut a >15cm piece from.
        l2_piece = l1_long_shrunk
        
        # Check if this piece can yield the 15cm target.
        # The target piece (15cm) must be the shorter piece from the second cut,
        # otherwise its shrunken length would not be an integer (15 * 0.75 = 11.25).
        if l2_piece > target_piece:
            l2_other_part = l2_piece - target_piece
            
            # The other part must be longer than the target piece,
            # and its length must be a multiple of 4 for its shrunken length to be an integer.
            if l2_other_part > target_piece and l2_other_part % 4 == 0:
                
                # Solution found. Print the detailed steps.
                print("Step-by-step derivation of the solution:")
                
                # First cut equation
                print(f"\n1. First, the {initial_length} cm rope is cut into a longer piece and a shorter piece.")
                print(f"   Equation: {l1_long} cm + {l1_short} cm = {initial_length} cm")

                # First shrink equation
                print(f"2. The longer piece ({l1_long} cm) shrinks by 25%.")
                print(f"   Equation: {l1_long} * {shrink_factor} = {l1_long_shrunk} cm")
                
                # Second cut equation
                print(f"3. Next, the {l1_long_shrunk} cm piece is cut to get our target piece.")
                print(f"   Equation: {l2_other_part} cm + {target_piece} cm = {l1_long_shrunk} cm")

                # Second shrink equation
                l2_other_shrunk = int(l2_other_part * shrink_factor)
                print(f"4. The longer portion of this second cut ({l2_other_part} cm) shrinks.")
                print(f"   Equation: {l2_other_part} * {shrink_factor} = {l2_other_shrunk} cm")

                print(f"\nThis proves it is possible to get a {target_piece} cm piece.")
                print(f"The length of the longer portion measured and cut from the initial rope was {l1_long} cm.")
                
                solution_found = True
                # Use a global or nonlocal variable if you want to pass this value out,
                # or just return it. For this script, we print and break.
                global final_answer
                final_answer = l1_long
                break

    if not solution_found:
        print("No solution was found that meets all the conditions.")

# Execute the function to find and print the solution.
solve_rope_puzzle()
<<<52>>>