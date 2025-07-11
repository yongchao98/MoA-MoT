def solve_rearrangement():
    """
    This function determines the substituents at the labeled positions
    after the backbone rearrangement and prints the result.
    """
    substituents = {
        1: 'CH3',
        2: 'H',
        3: 'CH3',
        4: 'CH3',
        5: 'H'
    }

    # Format the output string as "1 = V, 2 = W, 3 = X, 4 = Y, 5 = Z"
    output_parts = []
    for position in sorted(substituents.keys()):
        group = substituents[position]
        output_parts.append(f"{position} = {group}")

    print(", ".join(output_parts))

solve_rearrangement()