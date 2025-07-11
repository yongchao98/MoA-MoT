def solve_rearrangement():
    """
    Determines the substituents at the numbered positions in the product
    of the given rearrangement reaction.
    """
    substituents = {
        1: "CH3",
        2: "CH3",
        3: "H",
        4: "CH3",
        5: "H"
    }

    # Print the result in the specified format
    # The final equation is the list of substituents at each position.
    output_parts = []
    for position in sorted(substituents.keys()):
        substituent = substituents[position]
        output_parts.append(f"{position} = {substituent}")

    print(", ".join(output_parts))

solve_rearrangement()