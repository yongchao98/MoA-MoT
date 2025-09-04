import re

class OrganicChemistryProblemChecker:
    """
    A class to check the correctness of a multi-step organic chemistry problem.
    It validates the reasoning for each step of the synthesis and the final NMR analysis.
    """
    def __init__(self, question, llm_answer_text):
        self.question = question
        self.llm_answer_text = llm_answer_text.lower()
        self.final_answer_option = self.extract_final_answer(llm_answer_text)
        self.reasoning = self.extract_reasoning(llm_answer_text)
        self.options = {
            "A": "triplet",
            "B": "doublet of triplets",
            "C": "triplet of triplets",
            "D": "pentet"
        }

    def extract_final_answer(self, text):
        """Extracts the final answer option (e.g., 'A', 'B') from the text."""
        match = re.search(r'<<<([A-D])>>>', text)
        return match.group(1) if match else None

    def extract_reasoning(self, text):
        """Extracts the reasoning part of the answer."""
        return text.split('<<<')[0].strip().lower()

    def check_step_1_product_1(self):
        """
        Checks the identification of Product 1.
        Reaction: 1,3-dibromoadamantane + KOH, 240C -> Product 1
        Expected Product: Protoadamantan-4-one (via skeletal rearrangement)
        Constraints: IR at 1720 cm-1 (ketone) is consistent. The provided NMR data is a known red herring.
        """
        expected_product = "protoadamantan-4-one"
        if expected_product not in self.reasoning:
            return False, f"Incorrect identification of Product 1. The reaction of 1,3-dibromoadamantane with strong base at high temperature is a known rearrangement to form {expected_product}."
        return True, "Product 1 correctly identified."

    def check_step_2_product_2(self):
        """
        Checks the identification of Product 2.
        Reaction: Product 1 (ketone) + Al(O-iPr)3, heat -> Product 2
        Key Inference: The next step is ozonolysis, so Product 2 must be an alkene.
        Expected Logic: Meerwein-Ponndorf-Verley (MPV) reduction of the ketone to an alcohol, followed by heat-induced dehydration.
        Expected Product: Protoadamantene
        """
        expected_product = "protoadamantene"
        if expected_product not in self.reasoning:
            return False, f"Incorrect identification of Product 2. The reaction should lead to {expected_product} to allow for the subsequent ozonolysis step."
        
        # Check if the logic of reduction + dehydration is mentioned
        if not ("reduction" in self.reasoning and ("dehydration" in self.reasoning or "elimination" in self.reasoning)):
             return False, "Incorrect reasoning for Product 2 formation. The answer must mention both the reduction of the ketone and the subsequent dehydration to form the alkene (protoadamantene)."

        return True, "Product 2 correctly identified."

    def check_step_3_product_3(self):
        """
        Checks the identification of Product 3.
        Reaction: Product 2 (protoadamantene) + O3; then DMS -> Product 3
        Expected Logic: Reductive ozonolysis of the C=C bond in protoadamantene.
        Structure of Protoadamantene: The double bond is between two carbons each bearing one hydrogen.
        Expected Product: Bicyclo[3.3.1]nonane-3,7-dicarbaldehyde (a dialdehyde).
        """
        # Allow for variations like "dicarbaldehyde" or "dialdehyde"
        if "dicarbaldehyde" not in self.reasoning and "dialdehyde" not in self.reasoning:
            return False, "Incorrect functional group for Product 3. Ozonolysis of protoadamantene should yield a dialdehyde."
        
        if "bicyclo[3.3.1]nonane" not in self.reasoning:
            return False, "Incorrect skeleton for Product 3. Ozonolysis of protoadamantene yields a bicyclo[3.3.1]nonane skeleton."

        return True, "Product 3 correctly identified."

    def check_step_4_nmr_analysis(self):
        """
        Checks the NMR analysis of Product 3.
        1. Identify the most deshielded proton.
        2. Determine its coupling pattern.
        """
        # Check 1: Identification of the proton of interest.
        # The reasoning should acknowledge that aldehyde protons are most deshielded, but their pattern (doublet) is not an option,
        # thus shifting focus to the next most deshielded protons (H3/H7).
        if "aldehyde proton" not in self.reasoning:
            return False, "NMR analysis is incomplete. It fails to consider the aldehyde protons as the most deshielded."
        if not ("doublet" in self.reasoning and "not an option" in self.reasoning):
             return False, "NMR analysis is flawed. It should explain that the aldehyde proton's expected pattern (doublet) is not among the choices, which is a key step in the problem's logic."
        if "h3" not in self.reasoning and "h7" not in self.reasoning:
            return False, "NMR analysis is flawed. It fails to identify the next most deshielded protons (H3/H7) for analysis."

        # Check 2: Coupling pattern determination.
        # The logic should involve the dual-chair conformation, placing H3/H7 as axial,
        # and coupling to two axial and two equatorial neighbors.
        if "axial" not in self.reasoning or "equatorial" not in self.reasoning:
            return False, "NMR coupling analysis is incomplete. It must consider the axial/equatorial positions of neighboring protons in the chair conformation to determine the coupling constants."
        
        # The final pattern should be a triplet of triplets.
        expected_pattern = "triplet of triplets"
        if expected_pattern not in self.reasoning:
            return False, f"Incorrect coupling pattern deduced. The coupling of an axial proton (H3/H7) to two non-equivalent sets of two equivalent protons (axial and equatorial neighbors) results in a '{expected_pattern}'."
        
        return True, "NMR analysis is correct."

    def run_check(self):
        """
        Executes all checks and returns the final verdict.
        """
        if not self.final_answer_option:
            return "Could not find a final answer in the format <<<X>>>."

        correct_option = "C"
        
        # Step-by-step check of the reasoning
        checks = [
            self.check_step_1_product_1,
            self.check_step_2_product_2,
            self.check_step_3_product_3,
            self.check_step_4_nmr_analysis
        ]

        for check_func in checks:
            is_correct, reason = check_func()
            if not is_correct:
                return f"The reasoning is incorrect. Reason: {reason}"

        # If reasoning is sound, check if the final selected option matches the reasoning
        if self.final_answer_option != correct_option:
            correct_pattern = self.options[correct_option]
            return f"The reasoning correctly leads to '{correct_pattern}' (Option {correct_option}), but the selected answer was Option {self.final_answer_option}."

        return "Correct"

# The provided answer from the LLM to be checked
llm_answer_to_check = """
Based on a thorough analysis of the provided candidate answers and the underlying chemistry, here is a step-by-step breakdown to determine the correct solution.

**Step 1: Identify Product 1**
The reaction starts with 1,3-dibromoadamantane, which is heated with a strong base (KOH) at high temperature. The product data (IR: 1720 cm⁻¹ for a ketone; ¹H NMR: 14 protons total) points to a rearranged ketone isomer with the formula C₁₀H₁₄O. The most chemically sound identification, supported by known adamantane chemistry and proposed by multiple candidates (e.g., Answers 3, 4, 7, 8, 15), is **protoadamantan-4-one**. The unusual NMR data provided in the prompt is likely a red herring or erroneous, a common feature in such problems.

**Step 2: Identify Product 2**
Product 1 (protoadamantan-4-one) is heated with aluminum isopropoxide. The subsequent step is ozonolysis, which requires a carbon-carbon double bond. Therefore, the reaction with aluminum isopropoxide must be a two-step sequence occurring in one pot:
1.  **Meerwein-Ponndorf-Verley (MPV) Reduction:** The ketone is reduced to the corresponding alcohol, protoadamantan-4-ol.
2.  **Dehydration:** The alcohol eliminates water upon heating to form an alkene.
The resulting alkene, Product 2, is **protoadamantene**. This logical deduction is a key step and is made by the most credible candidate answers (e.g., 3, 4, 6, 7, 8).

**Step 3: Identify Product 3**
Product 2 (protoadamantene) undergoes reductive ozonolysis (O₃ followed by DMS). This reaction cleaves the C=C double bond. The double bond in protoadamantene is between two tertiary carbons (each bonded to one hydrogen). Therefore, the cleavage results in two aldehyde groups. The ring-opening of the protoadamantene skeleton yields **bicyclo[3.3.1]nonane-3,7-dicarbaldehyde**. This is the correct structure for Product 3, as identified in Answers 3 and 8. Other proposed structures like diketones or keto-aldehydes are incorrect based on the structure of protoadamantene.

**Step 4: Analyze the ¹H NMR Spectrum of Product 3**
The question asks for the coupling pattern of the most deshielded non-exchangeable proton in bicyclo[3.3.1]nonane-3,7-dicarbaldehyde.

1.  **Most Deshielded Protons:** The two aldehyde protons (-CHO) are the most deshielded (δ ≈ 9-10 ppm). Each would be coupled to the adjacent ring proton (H3 or H7), resulting in a doublet. Since "doublet" is not an option, we must consider the next most deshielded protons.
2.  **Next Most Deshielded Protons:** These are the protons on the carbons bearing the aldehyde groups, H3 and H7. They are alpha to a carbonyl group and are chemically equivalent due to the molecule's symmetry.
3.  **Coupling Analysis of H3/H7:** The bicyclo[3.3.1]nonane skeleton adopts a dual-chair conformation. To minimize steric hindrance, the bulky aldehyde groups occupy equatorial positions, forcing the H3 and H7 protons into axial positions. An axial proton like H3 is coupled to four neighboring protons on the adjacent methylene groups (C2 and C4):
    *   Two axial protons (H2ax, H4ax).
    *   Two equatorial protons (H2eq, H4eq).
4.  **Splitting Pattern:** Due to symmetry, the two axial neighbors are equivalent, and the two equatorial neighbors are equivalent. The coupling constant to the axial neighbors (J_ax) is large, while the coupling to the equatorial neighbors (J_eq) is small.
    *   The coupling to the two equivalent axial protons splits the signal into a **triplet**.
    *   Each peak of this triplet is then further split into another **triplet** by the two equivalent equatorial protons.
5.  **Conclusion:** The resulting pattern is a **triplet of triplets**.

This detailed analysis aligns perfectly with the reasoning in candidates 3 and 8, leading to the same conclusion.

<<<C>>>
"""

# Instantiate the checker and run the check
checker = OrganicChemistryProblemChecker(question="", llm_answer_text=llm_answer_to_check)
result = checker.run_check()
print(result)