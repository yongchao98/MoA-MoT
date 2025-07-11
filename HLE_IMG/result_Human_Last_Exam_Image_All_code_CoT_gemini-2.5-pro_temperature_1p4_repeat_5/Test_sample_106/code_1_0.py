def solve_rearrangement():
    """
    This function determines the substituents at positions 1, 2, 3, 4, and 5
    after the acid-catalyzed rearrangement of the given triterpenoid.
    """
    substituents = {
        1: "CH3",
        2: "H",
        3: "CH3",
        4: "CH3",
        5: "H"
    }

    # Print the result in the required format
    for position, substituent in substituents.items():
        print(f"{position} = {substituent}", end="")
        if position < 5:
            print(", ", end="")
    print()

solve_rearrangement()