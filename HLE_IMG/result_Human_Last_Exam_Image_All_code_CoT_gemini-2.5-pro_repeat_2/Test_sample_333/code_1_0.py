def solve_puzzle():
    """
    This function determines the correct pairing between the airfoil shape plots (A-H)
    and the pressure coefficient plots (1-8) based on aerodynamic principles.
    """
    
    # Based on the detailed analysis:
    # A (thin, moderate camber, rounded LE) -> 4 (moderate lift, rounded peak)
    # B (thin, moderate camber, sharp LE) -> 1 (moderate lift, sharp peak)
    # C (thin, low camber) -> 8 (low lift)
    # D (symmetrical) -> 2 (lowest lift)
    # E (thick, very blunt LE) -> 3 (moderate lift, very broad peak)
    # F (thick, high camber, rounded LE) -> 5 (high lift, broad peak)
    # G (thick, very high camber, sharp LE) -> 6 (highest lift, extreme peak)
    # H (thick, very high camber, sharp LE) -> 7 (very high lift, sharp peak)
    
    pairing = {
        'A': 4,
        'B': 1,
        'C': 8,
        'D': 2,
        'E': 3,
        'F': 5,
        'G': 6,
        'H': 7
    }
    
    # The final answer requires a sequence of eight integers without spaces or separators.
    # We get this by taking the values for keys A, B, C, D, E, F, G, H in order.
    
    answer_sequence = ""
    for airfoil_letter in sorted(pairing.keys()):
        answer_sequence += str(pairing[airfoil_letter])
        
    print(answer_sequence)

solve_puzzle()