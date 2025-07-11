def solve_puzzle():
    """
    This function determines the correct pairing between the airfoil shape plots (A-H)
    and the pressure coefficient plots (1-8).

    The reasoning is as follows:
    1.  Airfoils are categorized by thickness and camber, which determine lift potential.
    2.  Pressure plots are ranked by the enclosed area, which represents the generated lift.
    3.  The rankings are matched to find the correct pairings.

    -   A (thin, symmetric) -> 8 (lowest lift)
    -   B (thick, symmetric) -> 2 (low lift)
    -   C (thin, more camber) -> 1 (medium-high lift)
    -   D (thin, less camber) -> 5 (medium lift)
    -   E (thick, most camber) -> 6 (highest lift)
    -   F (thick, high camber) -> 7 (very high lift)
    -   G (thick, medium-high camber) -> 3 (high lift)
    -   H (thick, least camber in group) -> 4 (high lift)

    The sequence of plot numbers corresponding to airfoils A through H is derived from these pairings.
    """
    # The list 'pairings' stores the plot number for each airfoil from A to H.
    # A -> 8
    # B -> 2
    # C -> 1
    # D -> 5
    # E -> 6
    # F -> 7
    # G -> 3
    # H -> 4
    pairings = [8, 2, 1, 5, 6, 7, 3, 4]
    
    # The final answer is a sequence of these integers without separators.
    answer = "".join(map(str, pairings))
    print(answer)

solve_puzzle()