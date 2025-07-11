def calculate_belly_dance_turns():
    """
    Calculates the total number of right, left, and back turns in a
    FCBD® (FatChanceBellyDance®) sequence.
    """
    # Contributions from non-turning moves are 0.
    # We only need to sum the contributions from turning moves.

    # Swivel Step Half Turn (16 counts)
    # This involves two turns that result in facing the back.
    ssht_right = 0
    ssht_left = 0
    ssht_back = 2

    # Sunanda (once)
    # A full 360-degree turn presents the right side, back, and left side.
    sunanda_right = 1
    sunanda_left = 1
    sunanda_back = 1

    # Barrel Turn (once)
    # Another full 360-degree turn.
    barrel_right = 1
    barrel_left = 1
    barrel_back = 1

    # Summing the contributions for each direction
    # The other moves (Swivel Step, Balancing Step, Figure 8) contribute 0 turns.
    total_right = ssht_right + sunanda_right + barrel_right
    total_left = ssht_left + sunanda_left + barrel_left
    total_back = ssht_back + sunanda_back + barrel_back
    
    # Print the breakdown of the calculation
    print("Calculating the total turns for the sequence:")
    print(f"Total Right Turns = {ssht_right} (from Swivel Step Half Turn) + {sunanda_right} (from Sunanda) + {barrel_right} (from Barrel Turn) = {total_right}")
    print(f"Total Left Turns = {ssht_left} (from Swivel Step Half Turn) + {sunanda_left} (from Sunanda) + {barrel_left} (from Barrel Turn) = {total_left}")
    print(f"Total Back Turns = {ssht_back} (from Swivel Step Half Turn) + {sunanda_back} (from Sunanda) + {barrel_back} (from Barrel Turn) = {total_back}")
    
    # Print the final result in the required format
    print("\nFinal Count (Right/Left/Back):")
    print(f"{total_right}/{total_left}/{total_back}")


calculate_belly_dance_turns()