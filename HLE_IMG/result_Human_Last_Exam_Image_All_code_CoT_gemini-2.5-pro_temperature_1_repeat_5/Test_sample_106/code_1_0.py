def solve_substituents():
    """
    This function determines and prints the substituents at the numbered positions
    in the product molecule based on visual analysis of the provided image.
    """
    substituents = {
        1: "CH3",
        2: "CH3",
        3: "H",
        4: "CH3",
        5: "CH3"
    }

    # Format the output string as requested: "1 = V, 2 = W, 3 = X, 4 = Y, 5 = Z"
    output_parts = []
    for position in sorted(substituents.keys()):
        substituent = substituents[position]
        output_parts.append(f"{position} = {substituent}")
    
    final_answer_string = ", ".join(output_parts)
    print(final_answer_string)

solve_substituents()