def solve_rearrangement():
    """
    Determines the substituents at the numbered positions in the product of the
    Wagner-Meerwein rearrangement and prints the result.
    """
    substituents = {
        1: "CH3",
        2: "CH3",
        3: "H",
        4: "CH3",
        5: "CH3"
    }

    # Format the output string as requested
    output_parts = []
    for position in sorted(substituents.keys()):
        substituent = substituents[position]
        output_parts.append(f"{position} = {substituent}")

    final_answer = ", ".join(output_parts)
    print(final_answer)

solve_rearrangement()