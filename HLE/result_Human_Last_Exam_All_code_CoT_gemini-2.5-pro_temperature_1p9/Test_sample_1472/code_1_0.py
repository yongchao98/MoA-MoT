def determine_point_group():
    """
    This script explains the step-by-step determination of the symmetry point group
    for the molecule bis(2,5-dithiahexane)copper.
    """

    print("Step 1: Analyze the molecular structure.")
    print("  - The molecule is bis(2,5-dithiahexane)copper.")
    print("  - The central atom is Copper (Cu).")
    print("  - The ligand is 2,5-dithiahexane (CH3-S-CH2-CH2-S-CH3).")
    print("  - This is a bidentate ligand, forming a five-membered chelate ring (Cu-S-C-C-S).")
    print("  - The complex contains one Cu atom and two of these ligands, so the coordination number is 4.")
    print("-" * 40)

    print("Step 2: Determine the 3D geometry.")
    print("  - The complex, formally [Cu(C4H10S2)2]^n+, has a CuS4 core.")
    print("  - Experimental studies on this complex show that it adopts a distorted tetrahedral geometry.")
    print("  - The two chelate rings are not planar; they are puckered. In the most stable form, both rings have the same conformation.")
    print("-" * 40)

    print("Step 3: Identify the symmetry elements present.")
    print("  - E: The identity element is always present.")
    print("  - Cn (Proper Rotation Axis):")
    print("    - There is a C2 axis that passes through the copper atom and bisects the angle between the two ligands.")
    print("    - A 180-degree rotation around this axis exchanges the positions of the two identical ligands.")
    print("    - There are no other C2 axes perpendicular to this one, so it is not a D-type point group.")
    print("  - σ (Mirror Plane):")
    print("    - There is no horizontal mirror plane (σh) because the molecule is not planar.")
    print("    - There are no vertical (σv) or dihedral (σd) mirror planes. Any plane containing the C2 axis would cut through the ligands, but the ligands themselves are not symmetrical due to the puckered rings and methyl groups.")
    print("  - i (Inversion Center):")
    print("    - The molecule does not have a center of inversion. Inverting through the copper atom would not map the molecule onto itself.")
    print("  - Sn (Improper Rotation Axis):")
    print("    - No Sn axes are present other than S2 which is equivalent to i (which is absent).")
    print("-" * 40)

    print("Step 4: Assign the Point Group.")
    print("  - The only symmetry elements found are the identity (E) and a single C2 rotation axis.")
    print("  - A molecule that possesses only E and a Cn axis belongs to the Cn point group.")
    # The number in the point group is 2
    final_number = 2
    print(f"  - Since the only rotation axis is C{final_number}, the point group is C{final_number}.")
    print("-" * 40)


if __name__ == "__main__":
    determine_point_group()
