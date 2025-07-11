import sys

def find_point_group():
    """
    This script determines the symmetry point group for
    bis(2,5-dithiahexane)copper by following a logical, step-by-step process.
    """
    # 1. Define the molecule
    molecule_name = "bis(2,5-dithiahexane)copper"
    central_atom = "Cu (Copper)"
    ligand_name = "2,5-dithiahexane (dth)"
    ligand_formula = "CH3-S-CH2-CH2-S-CH3"

    # 2. Print initial information
    print(f"--- Determining the Point Group of {molecule_name} ---")
    print(f"Step 1: Analyze the Molecular Structure and Geometry")
    print("-----------------------------------------------------")
    print(f"The ligand, {ligand_name}, is bidentate, binding to the central {central_atom} atom via its two sulfur atoms.")
    print("'Bis' means there are two such ligands, giving a total coordination number of four.")
    print("A four-coordinate complex can be tetrahedral or square planar.")
    print("While Cu(II) complexes are often square planar, Cu(I) complexes are typically tetrahedral.")
    print("Without a specified oxidation state, we analyze the highest-symmetry, chemically plausible structure.")
    print("The crystal structure of the Cu(I) complex, [Cu(dth)2]+, shows a distorted tetrahedral geometry. We will analyze the idealized version of this structure.\n")

    # 3. Walk through the point group determination flowchart
    print("Step 2: Point Group Determination Flowchart")
    print("-------------------------------------------")

    # Q1: Linearity and High Symmetry
    print("-> Is the molecule linear or does it possess very high symmetry (Td, Oh, Ih)?")
    print("   No. The presence of two distinct, non-planar chelate rings reduces the symmetry.\n")

    # Q2: Principal axis (Cn)
    print("-> What is the highest order rotation axis (Cn)?")
    print("   Visualize the two ligands arranged around the central copper atom.")
    print("   An axis can be drawn through the copper atom that, upon a 180-degree rotation, perfectly interchanges the two ligands.")
    print("   This is a C2 axis. There are no axes of higher order (C3, C4, etc.).")
    print("   Therefore, the principal axis is a C2 axis, and n = 2.\n")

    # Q3: Perpendicular C2 axes
    print("-> Are there n=2 C2 axes perpendicular to the principal C2 axis?")
    print("   Yes. In the idealized tetrahedral arrangement, two additional C2 axes can be found.")
    print("   Each is perpendicular to the principal axis and to each other.")
    print("   A molecule with a Cn axis and n perpendicular C2 axes belongs to the D family of point groups.")
    print("   Thus, our molecule belongs to the D2 family.\n")

    # Q4: Mirror Planes
    print("-> Does the molecule have a horizontal mirror plane (σh)?")
    print("   No. A plane perpendicular to the principal axis would cut through the puckered ligands, so it cannot be a mirror plane.\n")

    print("-> Does the molecule have any vertical or dihedral mirror planes (σv or σd)?")
    print("   No. Due to the chirality (handedness) of the puckered chelate rings (C-S-C-C-S) and the attached methyl groups, no plane of symmetry can be found that contains the principal C2 axis.\n")

    # 4. Conclusion
    print("Step 3: Conclusion")
    print("------------------")
    print("The molecule has a C2 principal axis and two perpendicular C2 axes (D2 family).")
    print("It has no mirror planes of any kind (σh, σv, σd).")
    print("Therefore, the point group is D_n where n=2.")
    print("\nThe final point group for bis(2,5-dithiahexane)copper is composed of:")
    print("The character for the D family: D")
    print("The number representing the order of the principal axis: 2")

if __name__ == '__main__':
    find_point_group()