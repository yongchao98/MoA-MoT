def calculate_dance_turns():
    """
    Calculates the number of right, left, and back turns in a FCBD dance sequence.
    """
    # Initialize counters for each orientation.
    total_right_turns = 0
    total_left_turns = 0
    total_back_turns = 0

    # --- Analyze each move for turns ---

    # Swivel Step Half Turn (16 counts): Assumed to be a turn and return (full 360 rotation)
    # due to the long duration (standard is 8 counts for a half turn).
    # A full 360 rotation shows right, left, and back once.
    right_swivel_half_turn = 1
    left_swivel_half_turn = 1
    back_swivel_half_turn = 1
    
    # Sunanda (once): Two quick pivot turns. Shows right and left sides.
    # The 'back' is a transient position and not counted as a full presentation.
    right_sunanda = 1
    left_sunanda = 1
    back_sunanda = 0

    # Barrel Turn (once): A standard 360-degree turn.
    # Shows right, left, and back once.
    right_barrel_turn = 1
    left_barrel_turn = 1
    back_barrel_turn = 1
    
    # --- Sum the turns from all relevant moves ---
    
    # The other moves (Swivel Step, Balancing Step, Figure 8) do not involve rotation.
    # Their contribution is 0 for all turn types.
    
    total_right_turns = right_swivel_half_turn + right_sunanda + right_barrel_turn
    total_left_turns = left_swivel_half_turn + left_sunanda + left_barrel_turn
    total_back_turns = back_swivel_half_turn + back_sunanda + back_barrel_turn

    # --- Print the final breakdown and result ---
    print("Calculating total turns based on the dance moves:")
    
    # Print the equation for right side turns
    print(f"Right side turns = {right_swivel_half_turn} (Swivel Step Half Turn) + {right_sunanda} (Sunanda) + {right_barrel_turn} (Barrel Turn) = {total_right_turns}")

    # Print the equation for left side turns
    print(f"Left side turns = {left_swivel_half_turn} (Swivel Step Half Turn) + {left_sunanda} (Sunanda) + {left_barrel_turn} (Barrel Turn) = {total_left_turns}")

    # Print the equation for back turns
    print(f"Back turns = {back_swivel_half_turn} (Swivel Step Half Turn) + {back_sunanda} (Sunanda) + {back_barrel_turn} (Barrel Turn) = {total_back_turns}")
    
    print(f"\nThe final count is {total_right_turns} right side, {total_left_turns} left side, and {total_back_turns} back turns.")

calculate_dance_turns()