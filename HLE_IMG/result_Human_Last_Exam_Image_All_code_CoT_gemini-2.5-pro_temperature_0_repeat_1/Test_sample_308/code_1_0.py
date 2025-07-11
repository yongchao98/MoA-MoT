def solve_roller_puzzle():
    """
    This function determines the correct pairing between roller configurations (1-8)
    and their corresponding angular displacement plots (A-H).

    The logic is as follows:
    1. The number of lobes on the driving roller (green) corresponds to the number of
       periodic variations (wiggles) in the displacement plot.
    2. The magnitude of the speed variation (amplitude of the wiggles) depends on the
       extremeness of the roller shapes. Smoother shapes lead to smaller variations.

    - Plot A (2 cycles) -> Config 3 (2 lobes)
    - Plot B (6 cycles, large variation) -> Config 1 (6 lobes, most extreme shapes)
    - Plot C (4 cycles, smooth variation) -> Config 6 (4 lobes, smoothest shapes)
    - Plot D (3 cycles) -> Config 8 (3 lobes)
    - Plot E (5 cycles) -> Config 7 (5 lobes)
    - Plot F (6 cycles, small variation) -> Config 5 (6 lobes, smoothest shapes)
    - Plot G (4 cycles, complex variation) -> Config 4 (4 lobes, complex shapes)
    - Plot H (6 cycles, medium variation) -> Config 2 (6 lobes, intermediate shapes)

    The sequence is formed by listing the configuration number for each plot from A to H.
    """
    
    # The mapping from plot letter to configuration number
    pairings = {
        'A': 3,
        'B': 1,
        'C': 6,
        'D': 8,
        'E': 7,
        'F': 5,
        'G': 4,
        'H': 2
    }
    
    # Create the final sequence string by ordering the plot letters alphabetically
    final_sequence = ""
    for plot_letter in sorted(pairings.keys()):
        final_sequence += str(pairings[plot_letter])
        
    print(final_sequence)

solve_roller_puzzle()