def calculate_fcbd_turns():
    """
    Calculates the number of turns to the right, left, and back
    for a given FCBD® dance sequence.
    """
    print("Analyzing the FCBD® dance sequence step-by-step:\n")

    # Initialize contributions for each move
    # [right, left, back]
    swivel_step = [0, 0, 0]
    swivel_step_half_turn = [0, 0, 2] # 16 counts = 2 repetitions of an 8-count move, 1 back turn each
    sunanda = [0, 0, 1] # One turn that presents the back once
    balancing_step = [1, 1, 0] # Assumes one turn to the right side and one to the left side
    figure_8 = [0, 0, 0]
    barrel_turn = [1, 1, 1] # A 360 spin presents right, back, and left sides once

    # Store contributions in a list for the final equation
    contributions = {
        "Swivel Step": swivel_step,
        "Swivel Step Half Turn": swivel_step_half_turn,
        "Sunanda": sunanda,
        "Balancing Step": balancing_step,
        "Figure 8": figure_8,
        "Barrel Turn": barrel_turn
    }
    
    # Print the breakdown
    for move, turns in contributions.items():
        print(f"'{move}': {turns[0]} right turn(s), {turns[1]} left turn(s), {turns[2]} back turn(s).")
    
    print("\nCalculating the totals:\n")

    # Calculate and print the final equation for each direction
    right_nums = [c[0] for c in contributions.values()]
    left_nums = [c[1] for c in contributions.values()]
    back_nums = [c[2] for c in contributions.values()]

    total_right = sum(right_nums)
    total_left = sum(left_nums)
    total_back = sum(back_nums)
    
    right_eq = " + ".join(map(str, right_nums))
    left_eq = " + ".join(map(str, left_nums))
    back_eq = " + ".join(map(str, back_nums))

    print(f"Total Right Turns = {right_eq} = {total_right}")
    print(f"Total Left Turns = {left_eq} = {total_left}")
    print(f"Total Back Turns = {back_eq} = {total_back}")

    print(f"\nFinal Answer: The dancer turns her/his right side {total_right} times, left side {total_left} times, and back {total_back} times to the audience.")

calculate_fcbd_turns()