def solve_puzzle():
    """
    This function determines the correct pairing between the airfoil shape plots (A-H)
    and the pressure coefficient plots (1-8) based on aerodynamic principles.

    The matching is as follows:
    A (symmetrical, medium thickness) -> 1 (moderate lift, sharp peak)
    B (symmetrical, thick) -> 8 (moderate lift, broader peak)
    C (cambered, thin) -> 2 (high lift for a thin airfoil)
    D (symmetrical, very thin) -> 4 (lowest lift)
    E (very thick, high camber) -> 5 (high lift, broad suction peak)
    F (thinner, high camber) -> 6 (highest lift, extreme suction spike)
    G (thinner, medium camber) -> 7 (high lift, sharp peak)
    H (thinnest of group, low camber) -> 3 (lowest lift in the high-lift group)
    """
    
    # The sequence of plot numbers for airfoils A, B, C, D, E, F, G, H
    pairings = {
        'A': 1,
        'B': 8,
        'C': 2,
        'D': 4,
        'E': 5,
        'F': 6,
        'G': 7,
        'H': 3
    }
    
    # Create the final sequence string by ordering the plot numbers alphabetically by airfoil label
    final_sequence = ""
    for letter in sorted(pairings.keys()):
        final_sequence += str(pairings[letter])
        
    print(final_sequence)

solve_puzzle()