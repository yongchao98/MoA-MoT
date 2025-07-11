def solve_rearrangement():
    """
    Determines the substituents at the numbered positions after a Wagner-Meerwein rearrangement.
    The solution is based on the well-established mechanism for the acid-catalyzed backbone
    rearrangement of friedelinol to olean-12-ene.
    """

    # Dictionary to store the identified substituent for each position
    # The keys are the numbered positions in the product molecule.
    # The values are the substituents ('H' or 'CH3') that end up at these positions.
    substituents = {
        1: "CH3",  # The methyl group remaining at C-4
        2: "H",      # The hydride that shifted from C-9 to C-10
        3: "CH3",  # The methyl group that shifted from C-8 to C-9
        4: "CH3",  # The methyl group that shifted from C-14 to C-8
        5: "H"       # The hydride that shifted from C-13 to C-14
    }

    # Print the result in the specified format
    print(f"1 = {substituents[1]}, 2 = {substituents[2]}, 3 = {substituents[3]}, 4 = {substituents[4]}, 5 = {substituents[5]}")

solve_rearrangement()