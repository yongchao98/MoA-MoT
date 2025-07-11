def solve_rearrangement():
    """
    Determines the substituents at the specified positions based on the analysis
    of the acid-catalyzed rearrangement of a triterpenoid.
    """
    # The analysis concludes that the methyl groups on the main skeleton
    # are not involved in the rearrangement, based on the final product structure.
    # Therefore, the substituents correspond to those in a standard oleanane skeleton.
    substituents = {
        1: "CH3",  # Methyl group at C-4
        2: "CH3",  # Methyl group at C-10
        3: "H",    # Hydrogen atom at C-9
        4: "CH3",  # Methyl group at C-14
        5: "CH3"   # Methyl group at C-8
    }

    # Format the output as requested: 1 = V, 2 = W, 3 = X, 4 = Y, 5 = Z
    output_parts = []
    for position in sorted(substituents.keys()):
        substituent = substituents[position]
        output_parts.append(f"{position} = {substituent}")

    print(", ".join(output_parts))

solve_rearrangement()