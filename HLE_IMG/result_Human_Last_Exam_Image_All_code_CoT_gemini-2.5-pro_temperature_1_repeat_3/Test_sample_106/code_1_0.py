def solve_rearrangement():
    """
    This function determines the substituents at the numbered positions
    in the product of the given backbone rearrangement reaction.
    """

    # Based on the mechanistic analysis, we can assign the final substituent at each position.
    # The positions are mapped to their respective substituents (H or CH3).
    substituents = {
        1: "CH3",
        2: "H",
        3: "CH3",
        4: "CH3",
        5: "H"
    }

    # Prepare the list of strings for the final output.
    output_parts = []
    for position in sorted(substituents.keys()):
        substituent_identity = substituents[position]
        output_parts.append(f"{position} = {substituent_identity}")

    # Join the parts into a single string, separated by ", ".
    final_answer_string = ", ".join(output_parts)

    # Print the final result.
    print(final_answer_string)

solve_rearrangement()