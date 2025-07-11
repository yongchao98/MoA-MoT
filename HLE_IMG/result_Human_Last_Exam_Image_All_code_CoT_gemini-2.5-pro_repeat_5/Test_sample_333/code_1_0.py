def solve_airfoil_matching():
    """
    Solves the airfoil matching problem based on aerodynamic principles.

    The logic is as follows:
    1.  Airfoils are categorized by thickness (thin: A-D, thick: E-H) and camber.
    2.  Lift is proportional to the area of the Cp plots and increases with camber and thickness.
    3.  Match symmetric airfoils first (A, E), which have the lowest lift. The thicker one (E) has more lift than the thinner one (A).
        A (thin, sym) -> Plot 2 (lowest lift)
        E (thick, sym) -> Plot 8 (higher lift than 2)
    4.  Match the thin, cambered series (B < C < D in camber). Their lift should increase accordingly.
        B (low camber) -> Plot 1
        C (med camber) -> Plot 5
        D (high camber) -> Plot 4
        This matches the visual increase in area: Area(2) < Area(1) < Area(5) < Area(4).
    5.  Match the thick, cambered series (F < G < H in camber). These are the highest lift cases.
        F (low camber) -> Plot 3
        G (med camber) -> Plot 6
        H (high camber) -> Plot 7
        This matches the visual increase in area: Area(8) < Area(3) < Area(6) < Area(7).
    """

    # Pairings: Airfoil -> Cp Plot
    pairings = {
        'A': 2,
        'B': 1,
        'C': 5,
        'D': 4,
        'E': 8,
        'F': 3,
        'G': 6,
        'H': 7
    }

    # Generate the final answer string by ordering the plot numbers by airfoil letter (A-H)
    # The sorted() function on the dictionary keys ensures alphabetical order.
    answer_sequence = "".join([str(pairings[key]) for key in sorted(pairings.keys())])
    
    print(f"The final pairing from airfoil A to H is:")
    print(answer_sequence)

solve_airfoil_matching()