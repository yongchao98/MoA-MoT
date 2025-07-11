def solve_rearrangement():
    """
    Determines the substituents at the numbered positions after the rearrangement.
    This function encapsulates the reasoning explained above.
    """
    substituents = {
        1: "CH3",  # The remaining methyl group at C4.
        2: "H",    # A hydride from C9 migrates to C10.
        3: "CH3",  # A methyl group from C14 migrates to C8.
        4: "H",    # A hydride from C13 migrates to C14.
        5: "H"     # The hydrogen at C18, which is unaffected by the rearrangement.
    }

    # Print the result in the specified format
    output = []
    for position in sorted(substituents.keys()):
        output.append(f"{position} = {substituents[position]}")
    
    print(", ".join(output))

solve_rearrangement()