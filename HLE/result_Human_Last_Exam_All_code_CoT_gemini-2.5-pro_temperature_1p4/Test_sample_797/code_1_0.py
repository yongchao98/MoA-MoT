import collections

def solve_dance_turns():
    """
    Calculates the number of right, left, and back turns in a FCBD® sequence.
    """
    # Define the turn counts for each relevant move.
    # Assumes standard left turns.
    # Format: {"right": count, "left": count, "back": count}
    move_definitions = {
        "Swivel Step": {"right": 0, "left": 0, "back": 0},
        "Swivel Step Half Turn": {"right": 1, "left": 0, "back": 1},
        "Sunanda": {"right": 1, "left": 1, "back": 0},
        "Balancing Step": {"right": 0, "left": 0, "back": 0},
        "Figure 8": {"right": 0, "left": 0, "back": 0},
        "Barrel Turn": {"right": 1, "left": 1, "back": 0}
    }

    # The dance sequence with repetitions
    sequence = [
        ("Swivel Step", 1),
        ("Swivel Step Half Turn", 1),
        ("Sunanda", 1),
        ("Balancing Step", 1),
        ("Figure 8", 8),
        ("Barrel Turn", 1)
    ]

    totals = collections.defaultdict(int)
    components = collections.defaultdict(list)

    # Calculate total turns by summing up turns from each move in the sequence
    for move, repetitions in sequence:
        turns = move_definitions[move]
        for i in range(repetitions):
            if turns["right"] > 0:
                totals["right"] += turns["right"]
                components["right"].append(str(turns["right"]))
            if turns["left"] > 0:
                totals["left"] += turns["left"]
                components["left"].append(str(turns["left"]))
            if turns["back"] > 0:
                totals["back"] += turns["back"]
                components["back"].append(str(turns["back"]))

    print("To find the answer, we sum the turns for each rotational move in the sequence.")
    print("-----------------------------------------------------------------------------")

    # Build and print the equation for right side turns
    right_equation = " + ".join(components["right"])
    print(f"Right side turns come from: Swivel Step Half Turn, Sunanda, and Barrel Turn.")
    print(f"Calculation: {right_equation} = {totals['right']}")
    print("")

    # Build and print the equation for left side turns
    left_equation = " + ".join(components["left"])
    print(f"Left side turns come from: Sunanda and Barrel Turn.")
    print(f"Calculation: {left_equation} = {totals['left']}")
    print("")
    
    # Build and print the equation for back turns
    back_equation = " + ".join(components["back"])
    print(f"Back turns come from: Swivel Step Half Turn.")
    print(f"Calculation: {back_equation} = {totals['back']}")
    print("-----------------------------------------------------------------------------")

    print(f"Final Count (Right/Left/Back): {totals['right']}/{totals['left']}/{totals['back']}")

solve_dance_turns()
<<<C>>>