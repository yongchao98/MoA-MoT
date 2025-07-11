def solve_roller_puzzle():
    """
    Solves the roller drive configuration puzzle by matching configurations (1-8)
    to their displacement plots (A-H).
    """

    # The mapping is determined by analyzing the frequency and amplitude of variation
    # in the displacement plots and matching them to the geometry of the rollers.
    # Key: Plot Letter, Value: Configuration Number
    mapping = {
        'A': 3,
        'B': 6,
        'C': 8,
        'D': 2,
        'E': 5,
        'F': 4,
        'G': 1,
        'H': 7
    }

    # Justifications for each match based on analysis
    justifications = {
        'A': "Plot A has 2 large variation cycles, matching the 2 lobes of the driver in configuration 3.",
        'B': "Plot B has 4 smooth variation cycles, matching the 4-lobed driver and smooth elliptical driven roller in configuration 6.",
        'C': "Plot C has 3 distinct variation cycles, matching the 3-lobed driver in configuration 8.",
        'D': "Plot D has 6 sharp, spiky variation cycles, matching the 6-lobed driver and sharp triangular driven roller in configuration 2.",
        'E': "Plot E has 6 smooth variation cycles, matching the 6-lobed driver and the smooth 3-lobed driven roller in configuration 5.",
        'F': "Plot F has 8 small-amplitude variation cycles, matching the 8-featured driver and less extreme driven roller in configuration 4.",
        'G': "Plot G has 6 complex variation cycles, matching the 6-lobed driver and the moderately complex driven roller in configuration 1.",
        'H': "Plot H has 8 large-amplitude variation cycles, matching the 8-featured driver and sharp driven roller in configuration 7."
    }

    print("Matching roller configurations to displacement plots:")
    
    final_sequence = []
    # Iterate through plots A to H in alphabetical order
    for plot_letter in sorted(mapping.keys()):
        config_number = mapping[plot_letter]
        print(f"Plot {plot_letter} corresponds to configuration {config_number}. Reason: {justifications[plot_letter]}")
        final_sequence.append(str(config_number))

    # The final answer is the sequence of 8 integers without separators
    final_answer = "".join(final_sequence)
    print("\nThe final sequence of integers for plots A-H is:")
    print(final_answer)

solve_roller_puzzle()