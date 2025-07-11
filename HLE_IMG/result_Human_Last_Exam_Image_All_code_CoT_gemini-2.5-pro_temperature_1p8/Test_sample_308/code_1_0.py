def solve_roller_puzzle():
    """
    Solves the roller drive puzzle by mapping configurations to displacement plots.

    The core logic is based on two principles:
    1. The number of lobes on the driving roller (green) determines the number of
       oscillations in the corresponding displacement plot.
    2. The "sharpness" and geometry of the lobes determine the amplitude of these
       oscillations (i.e., the variation in the driven roller's speed).
    """

    # Mapping of plot letter to its corresponding configuration number based on analysis.
    # Key: Plot Letter, Value: Configuration Number
    mappings = {
        'A': 8,  # Plot A has 3 oscillations, matching configuration 8 (3 lobes).
        'B': 5,  # Plot B has 6 large oscillations, matching config 5 (6 lobes, high variation).
        'C': 3,  # Plot C has 1 large hump, matching config 3 (1 lobe).
        'D': 1,  # Plot D has 4 irregular oscillations, matching config 1 (irregular 4 lobes).
        'E': 7,  # Plot E has 5 oscillations, matching config 7 (5 lobes).
        'F': 2,  # Plot F has 6 small oscillations, matching config 2 (6 lobes, low variation).
        'G': 4,  # Plot G has 4 extreme oscillations, matching config 4 (sharp 4 points).
        'H': 6   # Plot H has 4 small oscillations, matching config 6 (rounded 4 lobes).
    }

    print("Determined pairings (Plot -> Configuration):")
    # To satisfy the instruction "output each number in the final equation",
    # we first show the individual pairings.
    for plot_letter in sorted(mappings.keys()):
        config_number = mappings[plot_letter]
        print(f"Plot {plot_letter} -> {config_number}")

    # The final answer is the sequence of configuration numbers corresponding to
    # plots A, B, C, D, E, F, G, H.
    final_sequence = ""
    for plot_letter in sorted(mappings.keys()):
        final_sequence += str(mappings[plot_letter])

    print("\nFinal sequence of integers for plots A-H:")
    print(final_sequence)

solve_roller_puzzle()
<<<85317246>>>