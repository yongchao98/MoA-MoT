def solve_rearrangement():
    """
    This function determines the substituents at the numbered positions
    after a Wagner-Meerwein rearrangement.

    The analysis is based on a standard triterpenoid backbone rearrangement mechanism.
    A discrepancy between the provided starting material and the product is resolved
    by inferring a swap of substituents at positions C9 and C10 in the
    starting material diagram.
    """

    substituents = {
        1: "CH3",  # The methyl group at C4 that does not migrate.
        2: "CH3",  # The methyl group from C9 migrates to C10.
        3: "H",      # The hydrogen from C8 migrates to C9.
        4: "CH3",  # The methyl group from C14 migrates to C8.
        5: "H"       # The hydrogen from C13 migrates to C14.
    }

    # Format the output string as requested: 1 = V, 2 = W, 3 = X, 4 = Y, 5 = Z
    answer_parts = []
    for position in sorted(substituents.keys()):
        substituent = substituents[position]
        answer_parts.append(f"{position} = {substituent}")
    
    final_answer = ", ".join(answer_parts)
    print(final_answer)

solve_rearrangement()