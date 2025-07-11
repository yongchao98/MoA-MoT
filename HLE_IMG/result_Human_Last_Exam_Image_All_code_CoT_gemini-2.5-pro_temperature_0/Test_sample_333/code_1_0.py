def solve_airfoil_matching():
    """
    This function determines the correct pairing between airfoil shapes (A-H)
    and their pressure coefficient distributions (1-8) based on aerodynamic principles.

    The logic is as follows:
    1.  Symmetric airfoils (D, H) produce the least lift and have characteristic Cp plots.
        - D is thin, causing a sharp suction peak -> Plot 2.
        - H is thick, causing a rounded suction peak -> Plot 8.
    2.  Thick, cambered airfoils (E, F, G) produce high lift with rounded suction peaks.
        Lift is proportional to camber (E > F > G).
        - E (high camber) -> Plot 3 (highest lift).
        - F (medium camber) -> Plot 7 (medium lift).
        - G (low camber) -> Plot 4 (lowest lift).
    3.  Thin, cambered airfoils (A, B, C) produce lift with sharp suction peaks.
        Lift is proportional to camber (A > B > C).
        - A (high camber) -> Plot 6 (highest lift).
        - B (medium camber) -> Plot 1 (medium lift).
        - C (low camber) -> Plot 5 (lowest lift).

    The final sequence is formed by listing the plot numbers corresponding to airfoils A through H.
    """
    
    # Pairings: Airfoil -> Cp Plot
    pairings = {
        'A': 6,
        'B': 1,
        'C': 5,
        'D': 2,
        'E': 3,
        'F': 7,
        'G': 4,
        'H': 8
    }
    
    # Generate the sequence string from A to H
    result_sequence = ""
    for airfoil_letter in sorted(pairings.keys()):
        result_sequence += str(pairings[airfoil_letter])
        
    print(f"The final sequence of plot numbers for airfoils A-H is:")
    print(result_sequence)
    print("\n<<<" + result_sequence + ">>>")

solve_airfoil_matching()