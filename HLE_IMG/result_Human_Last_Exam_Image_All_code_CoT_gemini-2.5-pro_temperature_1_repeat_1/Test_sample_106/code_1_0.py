def solve_rearrangement():
    """
    This function determines the substituents at positions 1-5 in the product
    of the Friedelane-Oleanane rearrangement.
    """
    
    # Based on tracing the [1,2] shifts in the backbone rearrangement:
    # Position 1 (at C4): Started with two CH3 groups. One migrated. One remains.
    # Position 2 (at C10): Started with one CH3. It migrated away. An H migrated in.
    # Position 3 (at C9): Started with one H. It migrated away. A CH3 migrated in.
    # Position 4 (at C8): Started with one CH3. It migrated away. Another CH3 migrated in.
    # Position 5 (at C14): Started with one CH3. It migrated away. An H migrated in.
    
    substituents = {
        1: "CH3",
        2: "H",
        3: "CH3",
        4: "CH3",
        5: "H"
    }

    # Format the output string as requested
    output_parts = []
    for position in sorted(substituents.keys()):
        substituent = substituents[position]
        output_parts.append(f"{position} = {substituent}")

    final_answer = ", ".join(output_parts)
    print(final_answer)

solve_rearrangement()