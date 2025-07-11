def calculate_fcbd_turns():
    """
    Calculates the number of turns in a given FCBD® dance sequence.
    """
    # Step 1: Define the number of turns for a single repetition of each move.
    # The interpretation is as follows:
    # - A "back turn" is counted when the turn's primary goal is to face the back.
    # - A Barrel Turn's pass-through of the back does not count.
    move_turns = {
        'Swivel Step': {'right': 0, 'left': 0, 'back': 0},
        'Swivel Step Half Turn': {'right': 1, 'left': 1, 'back': 1},
        'Sunanda': {'right': 1, 'left': 1, 'back': 1},
        'Balancing Step': {'right': 0, 'left': 0, 'back': 0},
        'Figure 8': {'right': 0, 'left': 0, 'back': 0},
        'Barrel Turn': {'right': 1, 'left': 1, 'back': 0}
    }

    # The sequence as described in the problem.
    # Note: Swivel Step (8 counts) is 1 rep.
    # Swivel Step Half Turn (16 counts) is 2 reps.
    sequence = [
        {'name': 'Swivel Step', 'reps': 1},
        {'name': 'Swivel Step Half Turn', 'reps': 2},
        {'name': 'Sunanda', 'reps': 1},
        {'name': 'Balancing Step', 'reps': 1},
        {'name': 'Figure 8', 'reps': 8},
        {'name': 'Barrel Turn', 'reps': 1}
    ]

    # Step 2: Calculate total turns and prepare equation strings
    totals = {'right': 0, 'left': 0, 'back': 0}
    equation_parts = {'right': [], 'left': [], 'back': []}

    for move in sequence:
        name = move['name']
        reps = move['reps']
        
        for direction in ['right', 'left', 'back']:
            turns_for_move = reps * move_turns[name][direction]
            if turns_for_move > 0:
                totals[direction] += turns_for_move
                # Add the number for this move to the equation string
                equation_parts[direction].append(str(turns_for_move))

    # Step 3: Format and print the output
    right_equation = " + ".join(equation_parts['right'])
    left_equation = " + ".join(equation_parts['left'])
    back_equation = " + ".join(equation_parts['back'])

    print(f"Calculating right side turns: {right_equation} = {totals['right']}")
    print(f"Calculating left side turns: {left_equation} = {totals['left']}")
    print(f"Calculating back turns: {back_equation} = {totals['back']}")
    print("-" * 20)
    print(f"The dancer turns their right side {totals['right']} times, left side {totals['left']} times, and back {totals['back']} times.")

calculate_fcbd_turns()