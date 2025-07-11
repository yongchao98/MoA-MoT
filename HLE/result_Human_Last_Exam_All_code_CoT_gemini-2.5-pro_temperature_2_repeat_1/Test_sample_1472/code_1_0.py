import textwrap

def get_point_group():
    """
    Determines and explains the point group of bis(2,5-dithiahexane)copper.
    """
    explanation = [
        "1. Molecular Composition:",
        "   - The central atom is Copper (Cu).",
        "   - The ligand is 2,5-dithiahexane (CH3-S-CH2-CH2-S-CH3). It is a bidentate ligand, meaning it binds to the copper atom at two points (the two sulfur atoms).",
        "   - The prefix 'bis' indicates that there are two of these ligands.",
        "\n2. Coordination and Geometry:",
        "   - Two bidentate ligands result in a four-coordinate complex: [Cu(MeSCH2CH2SMe)2].",
        "   - Four-coordinate complexes are typically either tetrahedral or square planar.",
        "   - For a d9 ion like Cu(II) with relatively bulky sulfur-donor ligands, the structure often distorts from a perfect geometry. Experimental studies on this and similar complexes show a distorted tetrahedral coordination around the copper atom.",
        "\n3. Symmetry Analysis:",
        "   - Let's model the complex as a generic tetrahedral [M(AA)2], where AA is a bidentate ligand.",
        "   - An idealized, highly symmetric tetrahedral structure that is merely squashed or elongated belongs to the D2d point group. This structure possesses a principal C2 axis, two perpendicular C2 axes, two dihedral mirror planes (σd), and an S4 improper rotation axis.",
        "   - However, in the real molecule, the two bidentate ligands wrap around the central copper atom, creating a chiral 'propeller' twist. The chelate rings (Cu-S-C-C-S) are also puckered.",
        "   - This inherent chirality and puckering remove the mirror planes (σd) and the S4 improper rotation axis from the D2d group.",
        "   - The symmetry elements that remain are:",
        "     - The identity element (E).",
        "     - Three perpendicular two-fold rotation axes (C2). One C2 axis relates the two separate ligands, while the other two C2 axes pass through the ligands.",
        "   - A point group containing only the identity element and three perpendicular C2 axes is the D2 point group.",
        "\n4. Conclusion:",
        "   - The combination of the tetrahedral coordination and the chiral arrangement of the two bidentate ligands reduces the symmetry from the parent D2d to D2."
    ]

    # Print the step-by-step explanation
    for line in explanation:
        print(textwrap.fill(line, width=80))

    # Print the final answer
    point_group = "D2"
    print(f"\nThe point group of bis(2,5-dithiahexane)copper is {point_group}.")

get_point_group()