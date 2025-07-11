def calculate_fcbd_turns():
    """
    Calculates the number of turns to the right side, left side, and back
    for a given FCBD® sequence.

    The logic is as follows:
    - A 1/4 turn (90 degrees) counts as 1 side turn.
    - A 1/2 turn (180 degrees) counts as 1 side turn and 1 back turn.
    - A full spin (360 degrees) is counted as showing each side once,
      but passing through the back is not counted as a separate back turn.
    """

    # --- Contributions from each move ---

    # Swivel Step (1/4 turn left): 1 right side
    swivel_step = {'right': 1, 'left': 0, 'back': 0}

    # Swivel Step Half Turn (1/2 turn left): 1 right side, 1 back
    swivel_half_turn = {'right': 1, 'left': 0, 'back': 1}

    # Sunanda (full 360 spin): 1 right side, 1 left side
    sunanda = {'right': 1, 'left': 1, 'back': 0}

    # Balancing Step: No turns
    balancing_step = {'right': 0, 'left': 0, 'back': 0}

    # Figure 8 (8 times): No turns
    figure_8 = {'right': 0, 'left': 0, 'back': 0}

    # Barrel Turn (full 360 spin): 1 right side, 1 left side
    barrel_turn = {'right': 1, 'left': 1, 'back': 0}

    # --- Summing the totals ---

    total_right = (swivel_step['right'] + swivel_half_turn['right'] + sunanda['right'] +
                   balancing_step['right'] + figure_8['right'] + barrel_turn['right'])

    total_left = (swivel_step['left'] + swivel_half_turn['left'] + sunanda['left'] +
                  balancing_step['left'] + figure_8['left'] + barrel_turn['left'])

    total_back = (swivel_step['back'] + swivel_half_turn['back'] + sunanda['back'] +
                  balancing_step['back'] + figure_8['back'] + barrel_turn['back'])

    # --- Printing the equations and final result ---
    
    print("Calculating the totals for right side/left side/back turns:\n")

    # Print Right Side Calculation
    print("Right Side Equation:")
    print(f"{swivel_step['right']} (from Swivel Step) + "
          f"{swivel_half_turn['right']} (from Swivel Step Half Turn) + "
          f"{sunanda['right']} (from Sunanda) + "
          f"{balancing_step['right']} (from Balancing Step) + "
          f"{figure_8['right']} (from Figure 8) + "
          f"{barrel_turn['right']} (from Barrel Turn) = {total_right}\n")

    # Print Left Side Calculation
    print("Left Side Equation:")
    print(f"{swivel_step['left']} (from Swivel Step) + "
          f"{swivel_half_turn['left']} (from Swivel Step Half Turn) + "
          f"{sunanda['left']} (from Sunanda) + "
          f"{balancing_step['left']} (from Balancing Step) + "
          f"{figure_8['left']} (from Figure 8) + "
          f"{barrel_turn['left']} (from Barrel Turn) = {total_left}\n")

    # Print Back Turn Calculation
    print("Back Turn Equation:")
    print(f"{swivel_step['back']} (from Swivel Step) + "
          f"{swivel_half_turn['back']} (from Swivel Step Half Turn) + "
          f"{sunanda['back']} (from Sunanda) + "
          f"{balancing_step['back']} (from Balancing Step) + "
          f"{figure_8['back']} (from Figure 8) + "
          f"{barrel_turn['back']} (from Barrel Turn) = {total_back}\n")

    print("-" * 30)
    print(f"Final Count: {total_right} right side / {total_left} left side / {total_back} back")
    print("-" * 30)

calculate_fcbd_turns()