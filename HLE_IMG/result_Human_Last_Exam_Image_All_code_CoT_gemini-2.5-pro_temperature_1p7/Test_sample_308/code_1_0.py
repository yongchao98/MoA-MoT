def solve_roller_puzzle():
    """
    Solves the roller drive matching puzzle.

    The function determines the correct pairing between roller configurations (1-8)
    and displacement plots (A-H) based on geometric analysis.

    The matching logic is as follows:
    1.  The number of lobes on the driver roller dictates the number of periodic
        variations (cycles) in the corresponding displacement plot.
    2.  The shape and amplitude of these variations are determined by the specific
        geometries of the driving and driven rollers. A higher radius ratio
        (driver/driven) at the contact point results in a steeper slope on the plot.

    Pairings based on lobe/cycle count:
    - Plot A (2 cycles) <-> Config 3 (2 lobes)
    - Plot B (5 cycles) <-> Config 2 (5 lobes)
    - Plot C (3 cycles) <-> Config 8 (3 lobes)
    - Plot E (6 cycles) <-> Config 5 (6 lobes)

    Pairings of 4-lobed configurations based on curve shape:
    - Plot F: Very smooth, regular cycles -> Config 6 (smoothest shapes)
    - Plot G: Complex cycles with a sub-feature -> Config 4 (complex driver shape)
    - Plot H: Most extreme variation in slope -> Config 1 (most extreme follower shape)
    - Plot D: Remaining 4-cycle plot -> Config 7 (remaining 4-lobe config)

    The final sequence is constructed by ordering the configuration numbers
    according to the alphabetical order of the plots (A through H).
    """

    # Dictionary mapping each plot (A-H) to its corresponding configuration (1-8)
    pairings = {
        'A': 3,
        'B': 2,
        'C': 8,
        'D': 7,
        'E': 5,
        'F': 6,
        'G': 4,
        'H': 1
    }

    # Generate the final answer string by concatenating the numbers in order from A to H
    final_sequence = ""
    # We iterate through the sorted keys ('A', 'B', 'C', ...) to ensure correct order
    for plot_letter in sorted(pairings.keys()):
        config_number = pairings[plot_letter]
        final_sequence += str(config_number)

    print(final_sequence)

solve_roller_puzzle()