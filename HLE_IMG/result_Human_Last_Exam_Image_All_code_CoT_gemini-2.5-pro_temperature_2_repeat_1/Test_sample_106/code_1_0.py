def solve_rearrangement():
    """
    This function determines the substituents at positions 1, 2, 3, 4, and 5
    after the acid-catalyzed rearrangement, based on the established mechanism for
    the friedelane-oleanane backbone rearrangement.
    """
    substituents = {
        1: "CH3",  # C-4 retains one methyl group.
        2: "H",    # C-10 methyl migrates and is replaced by an H from C-9.
        3: "CH3",  # C-8 H migrates and is replaced by the methyl from C-14.
        4: "CH3",  # C-14 ends up with a methyl group in the final oleanane skeleton.
        5: "H"     # The dashed position on C-14 is occupied by the H migrating from C-13.
    }

    # Print the answer in the specified format
    for position, substituent in substituents.items():
        print(f"{position} = {substituent}", end="")
        if position < 5:
            print(", ", end="")
    print()

solve_rearrangement()