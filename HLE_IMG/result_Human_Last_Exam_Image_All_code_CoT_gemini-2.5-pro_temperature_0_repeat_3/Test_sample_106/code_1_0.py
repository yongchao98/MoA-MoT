def solve_rearrangement():
    """
    This function determines the substituents at the numbered positions
    after the acid-catalyzed backbone rearrangement and prints the result.
    """
    # The substituents are determined by tracing the Wagner-Meerwein rearrangement mechanism.
    # 1: CH3 at C-4 is a spectator.
    # 2: H from C-9 migrates to C-10.
    # 3: CH3 from C-14 migrates to C-8.
    # 4: CH3 from C-8 migrates to C-9.
    # 5: H from C-13 migrates to C-14.
    substituents = {
        1: "CH3",
        2: "H",
        3: "CH3",
        4: "CH3",
        5: "H"
    }

    # Format the output string as requested.
    output_parts = []
    for position in sorted(substituents.keys()):
        substituent = substituents[position]
        output_parts.append(f"{position} = {substituent}")
    
    print(", ".join(output_parts))

solve_rearrangement()