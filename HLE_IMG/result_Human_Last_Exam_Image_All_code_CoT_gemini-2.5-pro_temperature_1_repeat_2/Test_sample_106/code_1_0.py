def solve_rearrangement():
    """
    This function determines the substituents at the numbered positions
    in the product of the given chemical reaction.

    The reaction is an acid-catalyzed dehydration followed by a cascade of
    1,2-Wagner-Meerwein shifts (a "backbone rearrangement").

    Mechanism breakdown and substituent tracking:
    1.  The -OH group leaves, forming a carbocation at C3.
    2.  A series of 1,2-shifts occur to stabilize the charge:
        - 1,2-CH3 shift from C4 to C3. Cation at C4.
        - 1,2-H shift from C5 to C4. Cation at C5.
        - 1,2-CH3 shift from C10 to C5. Cation at C10.
        - 1,2-H shift from C9 to C10. Cation at C9.
        - 1,2-CH3 shift from C8 to C9. Cation at C8.
        - 1,2-CH3 shift from C14 to C8. Cation at C14.
        - 1,2-CH3 shift from C13 to C14. Cation at C13.
    3.  A proton is eliminated from C12 to form the final double bond.

    Tracking the substituents to their final positions:
    - Position 1 (at C4): One of the two original CH3 groups from C4 remains.
    - Position 2 (at C10): The H from C9 migrates here.
    - Position 3 (at C5): The CH3 from C10 migrates here.
    - Position 4 (at C8): The CH3 from C14 migrates here.
    - Position 5 (at C9): The CH3 from C8 migrates here.
    """
    substituents = {
        1: "CH3",
        2: "H",
        3: "CH3",
        4: "CH3",
        5: "CH3"
    }

    # Format the output as requested: 1 = V, 2 = W, 3 = X, 4 = Y, 5 = Z
    output_parts = []
    for position in sorted(substituents.keys()):
        substituent = substituents[position]
        output_parts.append(f"{position} = {substituent}")

    print(", ".join(output_parts))

solve_rearrangement()