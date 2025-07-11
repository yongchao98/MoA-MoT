def solve_rope_puzzle():
    """
    This function solves the magical rope puzzle by programmatically searching for the correct sequence of cuts.
    """
    initial_length = 60
    target_length = 15
    shrink_factor = 0.75

    # l1 is the longer piece from the first cut. It must be > 30.
    # For l1 * shrink_factor to be an integer, l1 must be a multiple of 4.
    # So we iterate l1 from 32 up to 59, in steps of 4.
    for l1 in range(32, initial_length, 4):
        s1 = initial_length - l1
        # The longer piece shrinks.
        l1_shrunk = int(l1 * shrink_factor)
        
        # These are the two pieces after the first cut.
        pieces_after_first_cut = [s1, l1_shrunk]
        
        # Now, check if the second cut can be performed on either of these pieces.
        for L2 in pieces_after_first_cut:
            # To get the target_length, it must be the shorter piece from the second cut.
            # The other piece from this cut would be L2 - target_length.
            if L2 > target_length:
                longer_piece_second_cut = L2 - target_length
                
                # Check if the piece to be cut is longer than the target piece
                # and if the shrunken result will be an integer.
                if longer_piece_second_cut > target_length and longer_piece_second_cut % 4 == 0:
                    # Solution found.
                    print("Solution Found!\n")
                    print("Step-by-step breakdown:")
                    print(f"1. Initial rope length: {initial_length}cm.")
                    print(f"2. First cut is made, creating a longer piece of {l1}cm and a shorter piece of {s1}cm.")
                    print(f"3. The longer piece ({l1}cm) shrinks by 25%.")
                    print(f"   Equation: {l1} * {shrink_factor} = {l1_shrunk}")
                    print(f"4. The lengths of the two portions after the first cut are now {l1_shrunk}cm and {s1}cm.")
                    
                    longer_portion_after_first_cut = max(l1_shrunk, s1)
                    
                    print(f"\nThe longer of these two portions is {longer_portion_after_first_cut}cm.")
                    
                    # Explain the second cut for completeness
                    shorter_piece_second_cut = target_length
                    longer_piece_second_cut_shrunk = int(longer_piece_second_cut * shrink_factor)
                    print(f"\n5. To get the final {target_length}cm piece, the {L2}cm rope is cut.")
                    print(f"   It's cut into a longer {longer_piece_second_cut}cm piece and a shorter {shorter_piece_second_cut}cm piece.")
                    print(f"6. The {longer_piece_second_cut}cm piece shrinks to {longer_piece_second_cut_shrunk}cm, successfully yielding the desired {target_length}cm piece.")

                    print(f"\nFinal Answer: The length of the longer portion after the first cut is {longer_portion_after_first_cut}cm.")
                    return longer_portion_after_first_cut

    print("No solution found.")
    return None

# Execute the function to find and print the solution.
final_answer = solve_rope_puzzle()
# The final answer is requested in a specific format at the end of the response.
# The code above already prints the answer with explanation.
# This is just to fulfill the output format requirement.
if final_answer is not None:
    print(f"<<<{final_answer}>>>")
