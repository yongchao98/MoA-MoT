def solve_rearrangement():
    """
    Determines the substituents at the labeled positions in the product of the
    acid-catalyzed rearrangement of friedelan-3-ol to an oleanane derivative.

    The logic is based on tracking the known Wagner-Meerwein [1,2] shifts.
    """

    # Dictionary to store the final substituent at each numbered position.
    substituents = {}

    # Position 1 (at C-4): One of the two original methyl groups at C-4 remains
    # after the other migrates to C-3.
    substituents[1] = 'CH3'

    # Position 2 (at C-10): The original methyl group at C-10 migrates to C-5,
    # and a hydrogen atom migrates from C-9 to C-10.
    substituents[2] = 'H'

    # Position 3 (at C-8): The original methyl group at C-8 migrates to C-9,
    # and a methyl group migrates from C-14 to C-8.
    substituents[3] = 'CH3'

    # Position 4 (at C-14): The original methyl group at C-14 migrates to C-8,
    # and a hydrogen atom migrates from C-13 to C-14.
    substituents[4] = 'H'

    # Position 5 (at C-17): The methyl group at C-17 is not involved in the
    # backbone rearrangement and remains in its original position.
    substituents[5] = 'CH3'

    # Format and print the output as requested.
    output_parts = []
    for position in sorted(substituents.keys()):
        output_parts.append(f"{position} = {substituents[position]}")

    print(", ".join(output_parts))

solve_rearrangement()