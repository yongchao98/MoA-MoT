def solve_puzzle():
    """
    This function determines the correct pairing between the roller configurations and the displacement plots.
    
    The logic is as follows:
    1. The number of oscillations in a plot corresponds to the number of lobes on the driving (left, green) roller.
    2. The amplitude of the oscillations corresponds to how much the shapes deviate from a circle.
    
    Matching based on oscillation count:
    - Plot A (1 osc) -> Config 3 (1 lobe)
    - Plot D (3 osc) -> Config 8 (3 lobes)
    - Plot G (5 osc) -> Config 7 (5 lobes)
    
    Matching 4-oscillation plots (C, E, H) to 4-lobe configs (4, 6):
    - Config 4 has very spiky shapes, causing high amplitude oscillations -> Plot C.
    - Config 6 has smoother shapes, causing medium amplitude oscillations -> Plot E.
    
    Matching 6-oscillation plots (B, F) to 6-lobe configs (1, 2, 5):
    - Config 2's driven roller is pointier than Config 5's, causing higher amplitude oscillations -> Plot B.
    - Config 5 has a smoother driven roller -> Plot F.
    
    By elimination, the remaining pair must be Plot H and Config 1. This is the puzzle's trick, as Config 1 (6 lobes) produces Plot H (4 oscillations).
    
    The final sequence for plots A, B, C, D, E, F, G, H is derived from these pairings.
    """
    
    # Pairings: Plot -> Config
    pairings = {
        'A': 3,
        'B': 2,
        'C': 4,
        'D': 8,
        'E': 6,
        'F': 5,
        'G': 7,
        'H': 1
    }
    
    # Generate the sequence of numbers for plots A through H
    sequence = ""
    for plot_letter in sorted(pairings.keys()):
        sequence += str(pairings[plot_letter])
        
    print(sequence)

solve_puzzle()