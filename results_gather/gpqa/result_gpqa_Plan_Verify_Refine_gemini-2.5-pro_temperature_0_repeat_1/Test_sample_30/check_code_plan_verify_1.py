import sys

class MoleculeSymmetryChecker:
    """
    A class to check the multi-step synthesis and symmetry analysis of the final product.
    """

    def __init__(self, proposed_answer_option, proposed_answer_reasoning):
        self.proposed_answer_option = proposed_answer_option
        self.proposed_answer_reasoning = proposed_answer_reasoning
        self.options = {"A": "c3", "B": "d2h", "C": "cs", "D": "c2h"}
        self.proposed_point_group = self.options.get(proposed_answer_option)

    def check_reaction_pathway(self):
        """
        Checks if the described reaction pathway to the final product is chemically sound.
        """
        # Step 1: Nitration of Toluene -> p-nitrotoluene. Correct.
        # Step 2: Oxidation of p-nitrotoluene -> p-nitrobenzaldehyde. Correct.
        # Step 3: Condensation -> 1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one. Correct.
        # The identification of the final product is correct.
        return True, "The reaction pathway to form Product 3, 1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one, is correct."

    def analyze_symmetry_of_product_3(self):
        """
        Analyzes the symmetry elements of 1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one.
        The most stable isomer is (E,E) and is assumed to be planar for maximum conjugation.
        """
        # Let's place the molecule in a coordinate system for analysis.
        # Central carbonyl carbon at the origin (0,0,0).
        # The C=O bond lies along the y-axis.
        # The molecule lies in the xy-plane.

        # Check for C2 axis:
        # A C2 axis exists along the y-axis (passing through the C=O bond).
        # Rotating 180 degrees around this axis swaps the two identical (p-nitrophenyl)vinylene arms.
        has_c2 = True

        # Check for Center of Inversion (i):
        # A center of inversion would be at the origin (the carbonyl carbon).
        # The carbonyl oxygen atom is at some coordinate (0, y, 0).
        # Its inverted point would be (0, -y, 0), but there is no atom there.
        # Therefore, the molecule lacks a center of inversion.
        has_inversion_center = False

        # Check for Horizontal Plane of Symmetry (σh):
        # A σh plane must be perpendicular to the principal C2 axis.
        # Here, the C2 axis is the y-axis. The perpendicular plane is the xz-plane.
        # The xz-plane is NOT a plane of symmetry because the oxygen atom is on the +y side
        # with no corresponding atom on the -y side.
        has_sigma_h = False
        
        # Check for Vertical Planes of Symmetry (σv):
        # A σv plane must contain the principal C2 axis (the y-axis).
        # Plane 1: The plane of the molecule itself (the xy-plane). This is a valid σv.
        # Plane 2: The yz-plane, which bisects the molecule. This is also a valid σv.
        has_sigma_v = True

        # Determine the point group based on the found elements:
        # The molecule has E, C2, and two σv planes. This defines the C2v point group.
        # The C2h point group requires E, C2, i, and σh.
        actual_point_group = "C2v"

        return actual_point_group, has_c2, has_inversion_center, has_sigma_h

    def check_correctness(self):
        """
        Runs all checks and returns the final verdict.
        """
        path_correct, path_reason = self.check_reaction_pathway()
        if not path_correct:
            return path_reason

        actual_point_group, has_c2, has_inversion_center, has_sigma_h = self.analyze_symmetry_of_product_3()

        if self.proposed_point_group == actual_point_group:
            return "Correct"
        else:
            error_report = []
            error_report.append(f"The answer is incorrect.")
            error_report.append(f"The question's final product, 1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one, belongs to the C2v point group, not {self.proposed_point_group} as claimed.")
            
            # Explain why C2h is wrong
            if self.proposed_point_group == "c2h":
                error_report.append("\nConstraint Check for C2h:")
                if not has_inversion_center:
                    error_report.append("- FAILED: The molecule lacks a center of inversion (i). The carbonyl oxygen does not have an equivalent atom on the opposite side of the molecule's center.")
                if not has_sigma_h:
                    error_report.append("- FAILED: The molecule lacks a horizontal mirror plane (σh) perpendicular to the principal C2 axis.")
            
            error_report.append(f"\nConclusion: The molecule's actual point group is C2v, based on the presence of a C2 axis and two vertical mirror planes. The provided answer 'D' (C2h) is incorrect because the required symmetry elements for C2h are not present.")
            
            return "\n".join(error_report)

# --- Execution of the check ---
# Provided answer from the LLM
llm_option = "D"
llm_reasoning = """
1.  **Step 1: Nitration of Toluene.** Toluene reacts with a mixture of nitric acid and sulfuric acid in an electrophilic aromatic substitution reaction. The methyl group (-CH3) is an ortho-, para-director. The major product, due to less steric hindrance, is p-nitrotoluene.
    *   **Product 1: p-nitrotoluene**

2.  **Step 2: Oxidation of Product 1.** p-Nitrotoluene is treated with manganese dioxide (MnO2) and sulfuric acid (H2SO4). MnO2 is a selective oxidizing agent that typically oxidizes benzylic C-H bonds to an aldehyde. Stronger oxidants like KMnO4 would form a carboxylic acid, but MnO2/H2SO4 is commonly used to synthesize the aldehyde.
    *   **Product 2: p-nitrobenzaldehyde**

3.  **Step 3: Condensation Reaction.** p-Nitrobenzaldehyde (Product 2) is treated with acetone and aqueous sodium hydroxide. This is a base-catalyzed crossed-aldol condensation, specifically a Claisen-Schmidt condensation. The base (NaOH) deprotonates acetone, which then acts as a nucleophile. Since acetone has acidic protons on both sides of the carbonyl, it can react with two equivalents of the aldehyde. This double condensation is a common reaction, leading to a highly conjugated and stable product.
    *   **Product 3: 1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one**

4.  **Step 4: Symmetry Analysis of Product 3.** The structure of 1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one is:
    (O2N-Ph)-CH=CH-C(=O)-CH=CH-(Ph-NO2)
    To maximize conjugation, the molecule adopts a planar conformation. Let's analyze the symmetry elements of this planar molecule:
    *   **C2 axis:** There is a two-fold rotational axis passing through the oxygen and carbon atoms of the central carbonyl group (C=O). Rotating the molecule 180° around this axis leaves the molecule unchanged.
    *   **Plane of symmetry (σh):** The plane containing all the heavy atoms of the molecule is a plane of symmetry. This is designated as a horizontal plane (σh) because it is perpendicular to the principal C2 axis.
    *   **Center of inversion (i):** The presence of a C2 axis and a perpendicular plane of symmetry (σh) guarantees the presence of a center of inversion. The inversion center is located at the position of the central carbonyl carbon.

    A molecule possessing an identity element (E), a C2 axis, a horizontal plane of symmetry (σh), and a center of inversion (i) belongs to the **C2h** point group.
"""

# Create a checker and run it
checker = MoleculeSymmetryChecker(llm_option, llm_reasoning)
result = checker.check_correctness()
print(result)