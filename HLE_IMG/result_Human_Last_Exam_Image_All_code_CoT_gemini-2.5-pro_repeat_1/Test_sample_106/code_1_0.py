def solve_rearrangement():
    """
    Determines the substituents at the numbered positions after the
    acid-catalyzed backbone rearrangement of the triterpenoid alcohol.
    """

    # Based on the mechanism of the Wagner-Meerwein backbone rearrangement:
    substituents = {
        1: "CH3",  # The methyl group at C-4 that does not migrate.
        2: "H",      # The hydrogen that migrates from C-9 to C-10.
        3: "CH3",  # The methyl group that migrates from C-8 to C-9.
        4: "CH3",  # The methyl group that migrates from C-14 to C-8.
        5: "H"       # The hydrogen that migrates from C-13 to C-14.
    }

    # Format the output string as requested
    result = ", ".join([f"{pos} = {sub}" for pos, sub in substituents.items()])
    print(result)

solve_rearrangement()