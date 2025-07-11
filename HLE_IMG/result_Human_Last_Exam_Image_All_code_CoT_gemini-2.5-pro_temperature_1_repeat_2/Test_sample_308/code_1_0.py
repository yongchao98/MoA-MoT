def solve_roller_puzzle():
    """
    This function determines the pairing between roller configurations and displacement plots.
    Each integer in the final sequence represents the configuration plot number (1-8)
    that corresponds to each alphabetically ordered displacement plot (A-H).
    """

    # The logic for pairing is as follows:
    # The number of oscillations in a plot corresponds to the number of lobes on the driving (green) roller.
    # The amplitude of the oscillations corresponds to the magnitude of variation in the ratio of the radii (R_driver / R_driven).

    # Plot A: 3 complex oscillations -> Config 7 (3 complex lobes)
    # Plot B: 5 oscillations -> Config 2 (5 lobes)
    # Plot C: 2 oscillations -> Config 3 (2 lobes)
    # Plot D: 4 large oscillations -> Config 4 (4 lobes, large radius variation)
    # Plot E: 4 small oscillations -> Config 6 (4 lobes, small radius variation)
    # Plot F: 6 large oscillations -> Config 1 (6 lobes, large radius variation on driven roller)
    # Plot G: 6 small oscillations -> Config 5 (6 lobes, small radius variation on driven roller)
    # Plot H: 3 simple oscillations -> Config 8 (3 simple lobes)

    pairings = {
        'A': 7,
        'B': 2,
        'C': 3,
        'D': 4,
        'E': 6,
        'F': 1,
        'G': 5,
        'H': 8
    }

    # Generate the sequence of numbers for plots A through H.
    result_sequence = ""
    for plot_letter in sorted(pairings.keys()):
        result_sequence += str(pairings[plot_letter])

    print(result_sequence)

solve_roller_puzzle()