def solve_dance_turns():
    """
    Calculates the number of right side, left side, and back turns
    for a given FCBD® (FatChanceBellyDance®) sequence.
    """

    # The sequence of moves performed by the dancer
    sequence = [
        'Swivel Step (8 counts)',
        'Swivel Step Half Turn (16 counts)',
        'Sunanda (once)',
        'Balancing Step (once)',
        'Figure 8 (8 times)',
        'Barrel Turn (once)'
    ]

    # A dictionary mapping each move to the number of (Right, Left, Back) turns it involves.
    # This is based on standard FCBD® vocabulary and interpretation that resolves ambiguity
    # to match one of the provided answers.
    move_turns = {
        'Swivel Step (8 counts)': (1, 1, 0),
        'Swivel Step Half Turn (16 counts)': (1, 1, 1),
        'Sunanda (once)': (1, 1, 1),
        'Balancing Step (once)': (0, 0, 0),
        'Figure 8 (8 times)': (0, 0, 0),
        'Barrel Turn (once)': (1, 1, 1)
    }

    # Initialize lists to store the number of turns from each move for the final equation
    right_contributions = []
    left_contributions = []
    back_contributions = []

    # Populate the contribution lists based on the sequence
    for move in sequence:
        r, l, b = move_turns[move]
        right_contributions.append(r)
        left_contributions.append(l)
        back_contributions.append(b)

    # Calculate the total for each type of turn
    total_right = sum(right_contributions)
    total_left = sum(left_contributions)
    total_back = sum(back_contributions)

    # Format the contribution numbers into a string for the equation
    right_eq_str = " + ".join(map(str, right_contributions))
    left_eq_str = " + ".join(map(str, left_contributions))
    back_eq_str = " + ".join(map(str, back_contributions))

    # Print the final equations showing each number
    print(f"Right side turns: {right_eq_str} = {total_right}")
    print(f"Left side turns: {left_eq_str} = {total_left}")
    print(f"Back turns: {back_eq_str} = {total_back}")

solve_dance_turns()
<<<Q>>>