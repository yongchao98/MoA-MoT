def solve_dance_turns():
    """
    Calculates the total number of right, left, and back turns
    in the given FCBD® dance sequence.
    """
    # Initialize turn counts
    right_turns = 0
    left_turns = 0
    back_turns = 0
    
    # --- Analysis of each move ---

    # Swivel Step Half Turn (16 counts)
    # 16 counts implies the 8-count move is done twice.
    # Each time involves a turn to the back and a return to the front.
    # This results in 1 right, 1 left, and 1 back turn per sequence.
    swivel_half_turn_reps = 2
    swivel_half_turn_right = 1 * swivel_half_turn_reps
    swivel_half_turn_left = 1 * swivel_half_turn_reps
    swivel_half_turn_back = 1 * swivel_half_turn_reps

    right_turns += swivel_half_turn_right
    left_turns += swivel_half_turn_left
    back_turns += swivel_half_turn_back
    
    # Sunanda (once)
    # This move involves a turn to the back and a return to the front.
    sunanda_reps = 1
    sunanda_right = 1 * sunanda_reps
    sunanda_left = 1 * sunanda_reps
    sunanda_back = 1 * sunanda_reps

    right_turns += sunanda_right
    left_turns += sunanda_left
    back_turns += sunanda_back
    
    # Barrel Turn (once)
    # A full 360 spin. It presents both sides, but the back-facing
    # moment is transient and not counted as a 'presentation'.
    barrel_turn_reps = 1
    barrel_turn_right = 1 * barrel_turn_reps
    barrel_turn_left = 1 * barrel_turn_reps
    barrel_turn_back = 0 * barrel_turn_reps

    right_turns += barrel_turn_right
    left_turns += barrel_turn_left
    back_turns += barrel_turn_back
    
    # --- Final Calculation and Output ---
    
    print("Calculating the total number of turns (Right/Left/Back):")
    print("-" * 30)
    
    # Final equation for Right Turns
    print(f"Total Right Turns = {swivel_half_turn_right} (from Swivel Step Half Turn) + {sunanda_right} (from Sunanda) + {barrel_turn_right} (from Barrel Turn) = {right_turns}")
    
    # Final equation for Left Turns
    print(f"Total Left Turns = {swivel_half_turn_left} (from Swivel Step Half Turn) + {sunanda_left} (from Sunanda) + {barrel_turn_left} (from Barrel Turn) = {left_turns}")
    
    # Final equation for Back Turns
    print(f"Total Back Turns = {swivel_half_turn_back} (from Swivel Step Half Turn) + {sunanda_back} (from Sunanda) + {barrel_turn_back} (from Barrel Turn) = {back_turns}")
    
    print("-" * 30)
    print(f"The final counts are: {right_turns} right / {left_turns} left / {back_turns} back.")

solve_dance_turns()