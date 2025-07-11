def calculate_fcbd_turns():
    """
    Calculates the total number of turns in a given FCBD sequence.
    """
    # Initialize counts for each type of turn
    total_right_turns = 0
    total_left_turns = 0
    total_back_turns = 0

    # Define turns for each move in the sequence
    # Swivel Step Half Turn (16 counts) is interpreted as two 8-count Turkish Swivels (R and L lead)
    swivel_half_turn_r = 2
    swivel_half_turn_l = 2
    swivel_half_turn_b = 2

    # Sunanda (once) shows all three facings
    sunanda_r = 1
    sunanda_l = 1
    sunanda_b = 1

    # Barrel Turn (once) shows sides, but not the back due to the forward bend
    barrel_turn_r = 1
    barrel_turn_l = 1
    barrel_turn_b = 0
    
    # Other moves have 0 turns, so we only need to sum the turning moves.

    # Calculate total right turns
    total_right_turns = swivel_half_turn_r + sunanda_r + barrel_turn_r
    
    # Calculate total left turns
    total_left_turns = swivel_half_turn_l + sunanda_l + barrel_turn_l

    # Calculate total back turns
    total_back_turns = swivel_half_turn_b + sunanda_b + barrel_turn_b
    
    # Print the breakdown of the calculation
    print("Calculating the total turns for the sequence...")
    print(f"Right side turns = {swivel_half_turn_r} (from Swivel Step Half Turn) + {sunanda_r} (from Sunanda) + {barrel_turn_r} (from Barrel Turn) = {total_right_turns}")
    print(f"Left side turns  = {swivel_half_turn_l} (from Swivel Step Half Turn) + {sunanda_l} (from Sunanda) + {barrel_turn_l} (from Barrel Turn) = {total_left_turns}")
    print(f"Back turns       = {swivel_half_turn_b} (from Swivel Step Half Turn) + {sunanda_b} (from Sunanda) + {barrel_turn_b} (from Barrel Turn) = {total_back_turns}")
    
    print(f"\nFinal count (Right/Left/Back): {total_right_turns}/{total_left_turns}/{total_back_turns}")


calculate_fcbd_turns()