import collections

def solve_dance_turns():
    """
    Calculates the total number of right, left, and back turns in a
    FCBD® (FatChanceBellyDance®) sequence.
    """
    # --- Assumptions ---
    # 1. The dancer starts facing the audience (Front).
    # 2. A "turn" is counted each time a profile (right side, left side, back) 
    #    is presented to the audience during a rotation.
    # 3. The dancer's orientation is maintained between moves.
    # 4. For consistency, we assume turns are made over the right shoulder.

    # --- Initialization ---
    total_turns = collections.defaultdict(int)
    turn_contributions = collections.defaultdict(list)
    
    # --- Move-by-Move Analysis ---

    # Move 1: Swivel Step
    # No rotation. Orientation remains Front.
    
    # Move 2: Swivel Step Half Turn
    # 180-degree turn from Front. Path: Front -> Right Side -> Back.
    # Presents Right Side (1) and Back (1). Orientation becomes Back.
    total_turns['right'] += 1
    total_turns['back'] += 1
    turn_contributions['right'].append(1)
    turn_contributions['back'].append(1)

    # Move 3: Sunanda
    # 360-degree turn starting from Back. Path: Back -> Left Side -> Front -> Right Side -> Back.
    # Presents Left Side (1) and Right Side (1). Orientation remains Back.
    total_turns['right'] += 1
    total_turns['left'] += 1
    turn_contributions['right'].append(1)
    turn_contributions['left'].append(1)

    # Move 4: Balancing Step
    # No rotation. Orientation remains Back.

    # Move 5: Figure 8 (8 times)
    # No body rotation. Orientation remains Back.

    # Move 6: Barrel Turn
    # 360-degree turn starting from Back. Path: Back -> Left Side -> Front -> Right Side -> Back.
    # Presents Left Side (1) and Right Side (1). Orientation remains Back.
    total_turns['right'] += 1
    total_turns['left'] += 1
    turn_contributions['right'].append(1)
    turn_contributions['left'].append(1)

    # --- Final Calculation ---
    print("Calculating the total number of turns (Right/Left/Back) based on the sequence.")
    print("The moves contributing to the turns are:")
    print("- Swivel Step Half Turn: 1 Right, 1 Back")
    print("- Sunanda: 1 Right, 1 Left")
    print("- Barrel Turn: 1 Right, 1 Left")
    print("-" * 30)

    right_eq = " + ".join(map(str, turn_contributions['right']))
    left_eq = " + ".join(map(str, turn_contributions['left']))
    back_eq = " + ".join(map(str, turn_contributions['back']))

    print(f"Total Right Turns = {right_eq} = {total_turns['right']}")
    print(f"Total Left Turns  = {left_eq} = {total_turns['left']}")
    print(f"Total Back Turns  = {back_eq} = {total_turns['back']}")
    print("-" * 30)
    print(f"Final Count: {total_turns['right']} Right / {total_turns['left']} Left / {total_turns['back']} Back")

solve_dance_turns()
<<<C>>>