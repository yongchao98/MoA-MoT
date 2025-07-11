def solve_rearrangement():
    """
    Determines the substituents at the numbered positions after the rearrangement.
    """
    substituents = {
        1: "CH3",
        2: "H",
        3: "CH3",
        4: "H",
        5: "H"
    }

    # Print the result in the specified format
    output = []
    for position, substituent in sorted(substituents.items()):
        output.append(f"{position} = {substituent}")
    
    print(", ".join(output))

solve_rearrangement()