def solve_roller_puzzle():
    """
    This function determines and prints the sequence of roller configuration numbers
    that correspond to the displacement plots A through H.

    The matching is as follows:
    - Plot A: 3 large-amplitude wiggles -> Config 8
    - Plot B: 5 large-amplitude wiggles -> Config 4
    - Plot C: 3 small-amplitude wiggles -> Config 5
    - Plot D: 4 large-amplitude wiggles -> Config 6
    - Plot E: 5 small-amplitude wiggles -> Config 7
    - Plot F: 6 small-amplitude wiggles -> Config 2
    - Plot G: 6 large-amplitude wiggles -> Config 1
    - Plot H: 4 small-amplitude wiggles -> Config 3
    """

    # The sequence of configuration numbers (1-8) corresponding to plots A through H
    correspondence = {
        'A': 8,
        'B': 4,
        'C': 5,
        'D': 6,
        'E': 7,
        'F': 2,
        'G': 1,
        'H': 3
    }

    # Generate the final sequence string by ordering the plot labels alphabetically
    final_sequence = ""
    for plot_label in sorted(correspondence.keys()):
        final_sequence += str(correspondence[plot_label])

    print(final_sequence)

solve_roller_puzzle()