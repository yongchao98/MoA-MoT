def solve_rearrangement():
    """
    This function determines and prints the substituents at the numbered positions
    after the acid-catalyzed backbone rearrangement.
    """
    
    # Based on the analysis of the Wagner-Meerwein rearrangement mechanism:
    substituents = {
        1: "CH3", # At C4, one of the two methyls remains.
        2: "H",     # At C10, the methyl group leaves and a hydrogen moves in.
        3: "CH3", # At C5, the hydrogen leaves and a methyl group moves in.
        4: "CH3", # At C9, the hydrogen leaves and a methyl group moves in.
        5: "H"      # At C8, the methyl group leaves and a hydrogen moves in.
    }

    # Format the output as a single string: "1 = V, 2 = W, 3 = X, 4 = Y, 5 = Z"
    result_parts = []
    # Sort by key to ensure the order is 1, 2, 3, 4, 5
    for position in sorted(substituents.keys()):
        substituent = substituents[position]
        result_parts.append(f"{position} = {substituent}")
        
    final_answer_string = ", ".join(result_parts)
    print(final_answer_string)

solve_rearrangement()