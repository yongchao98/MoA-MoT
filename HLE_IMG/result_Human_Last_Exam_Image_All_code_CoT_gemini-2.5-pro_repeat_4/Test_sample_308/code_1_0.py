def solve_roller_drive():
    """
    This function determines the mapping between roller configurations and displacement plots
    and prints the resulting sequence.
    """

    # The determined pairings are:
    # A -> 3, B -> 1, C -> 5, D -> 4, E -> 2, F -> 7, G -> 6, H -> 8
    
    pairings = {
        'A': 3,
        'B': 1,
        'C': 5,
        'D': 4,
        'E': 2,
        'F': 7,
        'G': 6,
        'H': 8
    }

    print("The correspondence between displacement plots (A-H) and configuration plots (1-8) is as follows:")
    
    final_sequence = ""
    # Iterate through the plots in alphabetical order
    for plot_letter in sorted(pairings.keys()):
        config_number = pairings[plot_letter]
        # Output each number in the final equation
        print(f"Plot {plot_letter} corresponds to Configuration {config_number}")
        final_sequence += str(config_number)

    print("\nThe required sequence of eight integers is:")
    print(final_sequence)

solve_roller_drive()