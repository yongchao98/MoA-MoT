import re

class ChemistryLogicChecker:
    """
    This class encapsulates the logic for checking the multi-step organic chemistry problem.
    It uses a simplified knowledge base to verify the reaction sequence and then
    applies NMR splitting rules to check the final analysis.
    """
    def __init__(self):
        # A simplified knowledge base for this specific problem
        self.known_reactions = {
            "1,3-dibromoadamantane_KOH_heat": "protoadamantan-4-one",
            "protoadamantan-4-one_Al(OiPr)3_heat": "protoadamant-4-ene",
            "protoadamant-4-ene_ozonolysis": "cis-bicyclo[3.3.0]octane-3,7-dicarbaldehyde"
        }
        self.compound_properties = {
            "protoadamantan-4-one": {
                "functional_group": "ketone",
                "protons": 14
            },
            "protoadamant-4-ene": {
                "functional_group": "alkene"
            },
            "cis-bicyclo[3.3.0]octane-3,7-dicarbaldehyde": {
                "alpha_proton_neighbors": {
                    # For H3, neighbors are on C2 and C4. Each is a CH2 group.
                    "group_1_protons": 2,
                    "group_2_protons": 2
                }
            }
        }

    def check_step1(self, llm_reasoning):
        """Checks the identification of Product 1."""
        ir_peak = 1720
        total_protons = 14
        
        if "protoadamantan-4-one" not in llm_reasoning:
            return "Reasoning Error: The LLM failed to identify Product 1 as protoadamantan-4-one, which is the known product of this rearrangement."

        props = self.compound_properties["protoadamantan-4-one"]
        if not (1700 <= ir_peak <= 1725 and props["functional_group"] == "ketone"):
            return "Constraint Check Failed (Product 1): The IR peak at 1720 cm-1 indicates a ketone, but the proposed structure is not a ketone."
        if props["protons"] != total_protons:
            return f"Constraint Check Failed (Product 1): The 1H NMR shows {total_protons} protons, but the proposed structure has {props['protons']}."
        
        return None

    def check_step2(self, llm_reasoning):
        """Checks the identification of Product 2."""
        if "protoadamant-4-ene" not in llm_reasoning and "protoadamantene" not in llm_reasoning:
            return "Reasoning Error: The LLM failed to identify Product 2 as protoadamant-4-ene. The reaction is an MPV reduction followed by dehydration."
        
        props = self.compound_properties["protoadamant-4-ene"]
        if props["functional_group"] != "alkene":
            return "Constraint Check Failed (Product 2): The subsequent ozonolysis step requires Product 2 to be an alkene."
            
        return None

    def check_step3(self, llm_reasoning):
        """Checks the identification of Product 3."""
        if "cis-bicyclo[3.3.0]octane-3,7-dicarbaldehyde" not in llm_reasoning:
            return "Reasoning Error: The LLM failed to identify Product 3 as cis-bicyclo[3.3.0]octane-3,7-dicarbaldehyde, the correct ozonolysis product."
        return None

    def check_nmr_analysis(self, llm_reasoning, final_answer_char):
        """Checks the NMR coupling pattern analysis for Product 3."""
        # The most deshielded non-exchangeable protons are H3 and H7.
        # We analyze the neighbors of H3.
        neighbors = self.compound_properties["cis-bicyclo[3.3.0]octane-3,7-dicarbaldehyde"]["alpha_proton_neighbors"]
        
        # Correct reasoning:
        # H3 is coupled to two protons on C2 and two protons on C4.
        # Assuming conformational averaging, the two protons on C2 act as one group, splitting H3 into a triplet.
        # The two protons on C4 act as another group, splitting H3 into another triplet.
        # Since the environments of C2 and C4 are different, the coupling constants (J) are different.
        # Therefore, the pattern is a triplet of triplets.
        
        correct_pattern = "triplet of triplets"
        
        # Check if the LLM's reasoning aligns with this.
        if "triplet of triplets" not in llm_reasoning.lower():
             return f"Reasoning Error: The LLM's explanation for the splitting pattern is unclear or incorrect. The correct reasoning leads to a triplet of triplets."

        answer_map = {"A": "triplet of triplets", "B": "doublet of triplets", "C": "triplet", "D": "pentet"}
        llm_pattern = answer_map.get(final_answer_char)

        if llm_pattern != correct_pattern:
            return f"Incorrect Final Answer: The analysis points to a '{correct_pattern}', but the selected answer corresponds to '{llm_pattern}'."
            
        return None

    def check(self, llm_response):
        match = re.search(r'<<<([A-D])>>>', llm_response)
        if not match:
            return "Failure: Could not find the final answer in the standard format, e.g., <<<A>>>."
        
        final_answer_char = match.group(1)
        
        # Check the logical flow of the provided answer
        checks = [
            self.check_step1,
            self.check_step2,
            self.check_step3,
        ]
        
        for check_func in checks:
            error = check_func(llm_response)
            if error:
                return error
        
        # Final check of the NMR analysis against the answer
        error = self.check_nmr_analysis(llm_response, final_answer_char)
        if error:
            return error

        return "Correct"

# The user-provided LLM answer is passed to this function
def check_correctness(llm_answer):
    checker = ChemistryLogicChecker()
    result = checker.check(llm_answer)
    return result

# The user-provided LLM answer
llm_answer_string = """
This is a multi-step organic synthesis problem. I need to identify the structure of each product in the reaction sequence and then analyze the NMR spectrum of the final product.

**Plan:**
1.  Identify Product 1 from the reaction of 1,3-dibromoadamantane with KOH, using the provided spectral data (IR, NMR).
2.  Identify Product 2 from the reaction of Product 1 with aluminum isopropoxide. This reaction should produce an alkene, as the next step is ozonolysis.
3.  Identify Product 3 from the ozonolysis of Product 2.
4.  Analyze the structure of Product 3 to determine the coupling pattern of its most deshielded non-exchangeable proton.

**Step 1: Identify Product 1**
*   **Reactants:** 1,3-dibromoadamantane (C₁₀H₁₄Br₂) + excess KOH, 240°C.
*   **Product 1 Data:** IR at 1720 cm⁻¹ (ketone, C=O). ¹H NMR shows 14 protons total, meaning the formula is C₁₀H₁₄O.
*   **Analysis:** The reaction of 1,3-dihaloadamantanes under harsh basic conditions is a known skeletal rearrangement that produces a ketone. The specific product that matches the formula and spectral data is **protoadamantan-4-one** (tricyclo[4.3.1.0³,⁸]decan-4-one). The provided information about 1,3-dehydroadamantane describes a different possible reaction pathway, but the product of that pathway does not match the spectral data given in the problem. I will proceed based on the data in the prompt.
*   **Action:** Search to confirm the structure of the ozonolysis product of the alkene derived from protoadamantan-4-one, as this is the most complex step.

<search>ozonolysis of protoadamantene</search>

<information>
The reaction sequence starting from 1,3-dibromoadamantane is a known problem. The first product (Product 1) is indeed protoadamantan-4-one. The second reaction, MPV reduction with aluminum isopropoxide followed by dehydration, yields protoadamant-4-ene (Product 2). The subsequent reductive ozonolysis of protoadamant-4-ene (tricyclo[4.3.1.0^{3,8}]dec-4-ene) cleaves the C4=C5 double bond to form **cis-bicyclo[3.3.0]octane-3,7-dicarbaldehyde** (Product 3). The formula for this product is C₁₀H₁₄O₂, which is consistent with the starting material C₁₀H₁₄ (protoadamantene) undergoing ozonolysis.
</information>

**Plan Refinement:**
The search confirms the structures of the key intermediates and the final product.
1.  Product 1: Protoadamantan-4-one.
2.  Product 2: Protoadamant-4-ene.
3.  Product 3: cis-bicyclo[3.3.0]octane-3,7-dicarbaldehyde.
4.  Now, I will analyze the ¹H NMR spectrum of Product 3.

**Step 2: Analyze the ¹H NMR of Product 3**
*   **Structure:** Product 3 is cis-bicyclo[3.3.0]octane-3,7-dicarbaldehyde. This is a fused two-ring system made of two five-membered rings. The molecule has C₂ symmetry.
*   **Most Deshielded Proton:** The most deshielded protons are the aldehyde protons (-CHO), but they are typically excluded from such questions or would give a simple doublet. The next most deshielded non-exchangeable protons are those on the carbons attached to the aldehyde groups, H3 and H7. Due to symmetry, H3 and H7 are chemically equivalent.
*   **Coupling Analysis of H3:** We need to identify the neighbors of the H3 proton. H3 is on a CH group. It is adjacent to two CH₂ groups: C2 and C4.
    *   H3 is coupled to the two protons on C2.
    *   H3 is coupled to the two protons on C4.
*   **Splitting Pattern:** The five-membered rings in the bicyclo[3.3.0]octane system are flexible and undergo rapid conformational changes. This often averages the coupling constants of the two protons on an adjacent CH₂ group.
    *   The two protons on C2 will act as an equivalent pair, splitting the H3 signal into a triplet (n+1, where n=2). Let's call this coupling constant J_A.
    *   The two protons on C4 will also act as an equivalent pair, splitting the H3 signal into another triplet. Let's call this coupling constant J_B.
*   Since the chemical environments of C2 (next to a bridgehead) and C4 are different, the coupling constants will be different (J_A ≠ J_B).
*   Therefore, the signal for H3 is split into a triplet by the C2 protons, and each line of that triplet is further split into another triplet by the C4 protons. The resulting pattern is a **triplet of triplets**.

<<<A>>>
"""

# Run the check and print the result
result = check_correctness(llm_answer_string)
print(result)