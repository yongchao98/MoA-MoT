def calculate_fcbd_turns():
    """
    Calculates the number of right side, left side, and back turns
    in a given FCBD (FatChanceBellyDance) sequence.
    """

    # Turns for each move: [right, left, back]
    swivel_step = [0, 0, 0]
    # 16 counts = 2 repetitions of the move. Each rep involves a full rotation
    # and features the back being presented to the audience.
    swivel_step_half_turn = [2, 2, 2]
    # A 3/4 turn and 1/4 pivot back. A continuous turn that passes through
    # the back position but doesn't feature it.
    sunanda = [1, 1, 0]
    balancing_step = [0, 0, 0]
    figure_8 = [0, 0, 0]
    # A full 360 spin. A continuous turn that passes through the back position.
    barrel_turn = [1, 1, 0]

    # Summing the turns
    total_right = swivel_step[0] + swivel_step_half_turn[0] + sunanda[0] + balancing_step[0] + figure_8[0] + barrel_turn[0]
    total_left = swivel_step[1] + swivel_step_half_turn[1] + sunanda[1] + balancing_step[1] + figure_8[1] + barrel_turn[1]
    total_back = swivel_step[2] + swivel_step_half_turn[2] + sunanda[2] + balancing_step[2] + figure_8[2] + barrel_turn[2]

    # Printing the breakdown of the calculation
    print("Calculating the total number of turns for each orientation:")
    
    print("\nRight Side Turns:")
    print(f"{swivel_step_half_turn[0]} (from Swivel Step Half Turn) + {sunanda[0]} (from Sunanda) + {barrel_turn[0]} (from Barrel Turn) = {total_right}")

    print("\nLeft Side Turns:")
    print(f"{swivel_step_half_turn[1]} (from Swivel Step Half Turn) + {sunanda[1]} (from Sunanda) + {barrel_turn[1]} (from Barrel Turn) = {total_left}")

    print("\nBack Turns:")
    print(f"{swivel_step_half_turn[2]} (from Swivel Step Half Turn) + {sunanda[2]} (from Sunanda) + {barrel_turn[2]} (from Barrel Turn) = {total_back}")

    print(f"\nFinal Count (Right/Left/Back): {total_right}/{total_left}/{total_back}")

calculate_fcbd_turns()