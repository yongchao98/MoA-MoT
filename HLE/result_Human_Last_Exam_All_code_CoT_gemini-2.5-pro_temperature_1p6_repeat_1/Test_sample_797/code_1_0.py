def calculate_fcbd_turns():
    """
    Calculates the number of times a dancer shows their right side, left side,
    and back to the audience during a specific FCBD® sequence.
    
    The logic is:
    - A side (right/left) is counted if the dancer passes that orientation.
    - The back is only counted if a move *finishes* with the back to the audience.
    """

    # Define the turns for each move based on the logic.
    # Format: {"move_name": [Right_turns, Left_turns, Back_turns]}
    move_definitions = {
        "Swivel Step":                [0, 0, 0],
        "Swivel Step Half Turn":      [1, 0, 1],
        "Sunanda":                    [1, 1, 0],
        "Balancing Step":             [0, 0, 0],
        "Figure 8":                   [0, 0, 0],
        "Barrel Turn":                [1, 1, 0],
    }

    # The sequence performed by the dancer
    sequence = [
        ("Swivel Step", 1),
        ("Swivel Step Half Turn", 1),
        ("Sunanda", 1),
        ("Balancing Step", 1),
        ("Figure 8", 8),
        ("Barrel Turn", 1),
    ]

    total_turns = {"Right": 0, "Left": 0, "Back": 0}
    
    # Lists to store the numbers for the final equation string
    right_components = []
    left_components = []
    back_components = []

    for move, repetitions in sequence:
        turns = move_definitions[move]
        
        # We only care about non-zero turn contributions for the equation
        if turns[0] > 0:
            right_components.append(str(turns[0] * repetitions))
        if turns[1] > 0:
            left_components.append(str(turns[1] * repetitions))
        if turns[2] > 0:
            back_components.append(str(turns[2] * repetitions))

        total_turns["Right"] += turns[0] * repetitions
        total_turns["Left"] += turns[1] * repetitions
        total_turns["Back"] += turns[2] * repetitions

    # Create the equation strings
    right_eq = " + ".join(right_components) if right_components else "0"
    left_eq = " + ".join(left_components) if left_components else "0"
    back_eq = " + ".join(back_components) if back_components else "0"

    print("Calculation of turns:")
    print(f"Right side to audience: {right_eq} = {total_turns['Right']}")
    print(f"Left side to audience: {left_eq} = {total_turns['Left']}")
    print(f"Back to audience: {back_eq} = {total_turns['Back']}")
    
    print(f"\nThe final count is {total_turns['Right']}/{total_turns['Left']}/{total_turns['Back']}.")

calculate_fcbd_turns()