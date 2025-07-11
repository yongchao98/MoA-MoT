def solve_rearrangement():
    """
    Determines the substituents at labeled positions after the acid-catalyzed
    rearrangement of β-amyrin to olean-12-ene.
    
    The reaction proceeds via a Wagner-Meerwein backbone rearrangement, which
    involves a cascade of 1,2-hydride and 1,2-methyl shifts. By tracing the
    movement of these groups, we can identify the final substituent at each
    labeled position.
    """

    # Analysis of the rearrangement cascade:
    # - H from C9 moves to C10 (Position 1)
    # - CH3 from C8 moves to C9 (Position 2)
    # - CH3 from C14 moves to C8 (Position 3)
    # - H from C13 moves to C14 (Position 4)
    # - C13 becomes part of the C=C bond. The H at C18 (near Position 5) is unaffected.

    substituents = {
        1: "H",
        2: "CH3",
        3: "CH3",
        4: "H",
        5: "H"
    }

    # Format the output string as requested
    output = []
    for position, substituent in substituents.items():
        output.append(f"{position} = {substituent}")
    
    print(", ".join(output))

solve_rearrangement()