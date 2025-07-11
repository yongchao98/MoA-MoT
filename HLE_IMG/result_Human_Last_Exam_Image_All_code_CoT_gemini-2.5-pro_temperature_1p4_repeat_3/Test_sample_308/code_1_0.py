def solve_roller_puzzle():
    """
    This function determines and prints the matching sequence for the roller drive puzzle.
    """
    # The dictionary stores the final determined pairings.
    # Key: Displacement plot letter (A-H)
    # Value: Roller configuration number (1-8)
    pairings = {
        'A': 1,
        'B': 2,
        'C': 3,
        'D': 4,
        'E': 8,
        'F': 5,
        'G': 7,
        'H': 6,
    }

    print("Roller Drive Configuration and Displacement Plot Pairings:")
    print("---------------------------------------------------------")
    final_sequence = []
    # Loop through the plots in alphabetical order and print the corresponding configuration.
    for plot_letter in sorted(pairings.keys()):
        config_number = pairings[plot_letter]
        print(f"Displacement Plot {plot_letter} corresponds to Configuration {config_number}")
        final_sequence.append(str(config_number))

    # The final answer is the sequence of configuration numbers for plots A through H.
    final_answer = "".join(final_sequence)
    print("\nFinal sequence of configuration numbers for plots A-H:")
    print(final_answer)

solve_roller_puzzle()