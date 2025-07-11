def solve_scansion():
    """
    This function scans the Latin line "et tibi bene esse soli quom sibi sit male".

    The scansion is based on the iambic senarius meter, which requires assuming
    contractions common in early Roman dramatic poetry to fit the 12-syllable pattern.
    L = Long syllable
    S = Short syllable
    """

    # Define the six metrical feet of the iambic senarius line.
    # Each foot is a tuple containing the syllables and their determined length.
    foot1 = ("et", "ti", "S", "S")
    foot2 = ("bi", "ben(e)", "S", "L")
    foot3 = ("(es)se", "sol(i)", "S", "L")
    foot4 = ("quom", "si(bi)", "L", "S")
    foot5 = ("si)bi", "sit", "S", "L")
    foot6 = ("ma", "le", "S", "L") # Final 'L' by Brevis in Longo

    feet = [foot1, foot2, foot3, foot4, foot5, foot6]

    final_scansion_parts = []
    for foot in feet:
        # We need to output each character of the scansion individually,
        # followed by a space for the foot boundary.
        final_scansion_parts.append(foot[2])
        final_scansion_parts.append(" ")
        final_scansion_parts.append(foot[3])
        final_scansion_parts.append("  ") # Space between feet

    # The prompt requires printing each symbol in the final equation.
    # We will print the letters for the scansion one by one.
    print_list = []
    for part in final_scansion_parts:
        for char in part:
            print_list.append(char)
            
    print("".join(print_list).strip())


solve_scansion()