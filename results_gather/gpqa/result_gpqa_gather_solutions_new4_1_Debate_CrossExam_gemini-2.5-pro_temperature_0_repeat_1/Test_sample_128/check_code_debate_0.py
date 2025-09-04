import re

class ChemistryProblemChecker:
    """
    This class checks the correctness of the answer to a multi-step organic chemistry problem.
    It simulates the reaction sequence logically to verify the final product.
    """

    def __init__(self, question_details, llm_answer):
        """
        Initializes the checker with the problem's constraints and the provided answer.
        
        Args:
            question_details (dict): A dictionary containing the options and hints.
            llm_answer (str): The full text of the answer provided by the LLM.
        """
        self.options = question_details.get("options", {})
        self.hints = question_details.get("hints", {})
        self.llm_answer_text = llm_answer

    def get_ring_size_from_name(self, name):
        """Determines the ring size of a cyclic compound from its IUPAC name."""
        if "cyclobutan" in name: return 4
        if "cyclopentan" in name: return 5
        if "cyclohexan" in name: return 6
        if "cycloheptan" in name: return 7
        return None

    def check_ir_consistency(self, compound_name, ir_value):
        """
        Checks if the IR value from the hint is consistent with the compound's ring size.
        - Strained 5-membered ring ketones absorb at higher frequencies (~1750 cm^-1).
        - Less strained 6-membered ring ketones absorb at lower frequencies (~1715 cm^-1).
        """
        ring_size = self.get_ring_size_from_name(compound_name)
        if ring_size == 5 and ir_value == 1750:
            return True, ""
        if ring_size == 6 and ir_value == 1715:
            return True, ""
        
        expected_ir = "unknown"
        if ring_size == 5: expected_ir = "~1750 cm^-1"
        if ring_size == 6: expected_ir = "~1715 cm^-1"
        
        return False, f"IR mismatch for {compound_name}. Expected IR for a {ring_size}-membered ring ketone is {expected_ir}, but hint gives {ir_value} cm^-1."

    def identify_compound_A(self):
        """
        Identifies Compound A using Hint (a) and verifies with Hint (b).
        - Hint (a) describes a Wittig reaction. A retro-Wittig analysis on the product
          `1,2-dimethyl-4-(propan-2-ylidene)cyclopentane` reveals the starting ketone.
        - The product implies a ketone at position 4 of a 1,2-dimethylcyclopentane.
        - IUPAC rules require numbering the carbonyl as C1, which places the methyls at C3 and C4.
        """
        compound_A_name = "3,4-dimethylcyclopentan-1-one"
        
        # Verify with Hint (b): IR of A is ~1750 cm^-1
        is_consistent, reason = self.check_ir_consistency(compound_A_name, 1750)
        if not is_consistent:
            return None, f"Hint (b) for Compound A is inconsistent with the structure derived from Hint (a). {reason}"
            
        return compound_A_name, ""

    def perform_tiffeneau_demjanov(self, starting_ketone):
        """
        Simulates the Tiffeneau-Demjanov rearrangement.
        This sequence (ketone -> cyanohydrin -> amino alcohol -> rearrangement)
        results in a one-carbon ring expansion of a cycloalkanone.
        """
        if "cyclo" not in starting_ketone or "one" not in starting_ketone:
            return None, f"Starting material '{starting_ketone}' is not a cycloalkanone, so Tiffeneau-Demjanov rearrangement is not the expected pathway."

        ring_size = self.get_ring_size_from_name(starting_ketone)
        if ring_size is None:
            return None, f"Could not determine ring size for '{starting_ketone}'."

        new_ring_size = ring_size + 1
        ring_map = {5: "cyclopentan", 6: "cyclohexan", 7: "cycloheptan"}
        
        if new_ring_size not in ring_map:
            return None, f"Ring expansion resulted in an unsupported ring size: {new_ring_size}"
            
        new_ring_name = ring_map[new_ring_size]

        # The substituent pattern (3,4-dimethyl) is maintained relative to the carbonyl.
        final_product_name = f"3,4-dimethyl{new_ring_name}-1-one"
        
        return final_product_name, ""

    def check(self):
        """
        Executes the full checking process.
        1. Extracts the chosen answer from the LLM's response.
        2. Deduces the correct product by simulating the reaction sequence.
        3. Verifies all steps against the provided hints.
        4. Compares the deduced correct answer with the LLM's chosen answer.
        """
        # Step 1: Extract the final answer choice from the LLM's response
        match = re.search(r'<<<([A-D])>>>', self.llm_answer_text)
        if not match:
            return "Incorrect: The final answer is not in the required format '<<<A>>>', '<<<B>>>', '<<<C>>>', or '<<<D>>>'."
        
        chosen_option = match.group(1)
        chosen_answer_name = self.options.get(chosen_option)

        # Step 2: Deduce the structure of Compound A from the hints
        compound_A, reason_A = self.identify_compound_A()
        if not compound_A:
            return f"Incorrect: {reason_A}"

        # Step 3: Simulate the reaction sequence to find Compound E
        compound_E, reason_E = self.perform_tiffeneau_demjanov(compound_A)
        if not compound_E:
            return f"Incorrect: {reason_E}"

        # Step 4: Verify the deduced Compound E with Hint (b)
        is_consistent, reason_ir_E = self.check_ir_consistency(compound_E, 1715)
        if not is_consistent:
            return f"Incorrect: The derived structure for Compound E ({compound_E}) is inconsistent with Hint (b). {reason_ir_E}"

        # Step 5: Compare the deduced correct answer with the LLM's chosen answer
        if compound_E == chosen_answer_name:
            return "Correct"
        else:
            correct_option_list = [key for key, val in self.options.items() if val == compound_E]
            if not correct_option_list:
                 return f"Incorrect: The deduced correct product is '{compound_E}', which is not among the options."
            
            correct_option = correct_option_list[0]
            return (f"Incorrect: The chemical reasoning correctly leads to Compound E being '{compound_E}' (Option {correct_option}). "
                    f"However, the final answer provided was Option {chosen_option} ('{chosen_answer_name}').")

# --- Execution ---
# Define the problem's parameters based on the provided text
problem_details = {
    "options": {
        "A": "2,2,3,4-tetramethylcyclobutan-1-one",
        "B": "3,4-dimethylcyclohexan-1-one",
        "C": "4-methylcycloheptan-1-one",
        "D": "2,3,4-trimethylcyclopentan-1-one"
    },
    "hints": {
        "a": "Compound A, when reacted with an unknown phosphorous ylide, forms 1,2-dimethyl-4-(propan-2-ylidene)cyclopentane.",
        "b": "The IR spectrum of compound A exhibits a strong peak at around 1750 cm^-1 whereas, IR spectrum of compound E exhibits strong peak at around 1715 cm^-1."
    }
}

# This is the final answer provided by the LLM that needs to be checked.
llm_final_answer = """
Here is a step-by-step analysis to determine the final answer.

### **Analysis of the Reaction Scheme**

The core of this problem is to correctly identify the starting material, Compound A, and then follow a known reaction pathway, the Tiffeneau-Demjanov rearrangement, to find the final product, Compound E.

**1. Identifying Compound A:**
*   **Hint (a)** describes a Wittig reaction. In this reaction, a ketone's carbonyl group (C=O) is replaced by a carbon-carbon double bond from a phosphorus ylide.
*   The product is `1,2-dimethyl-4-(propan-2-ylidene)cyclopentane`. The `(propan-2-ylidene)` group is `=C(CH₃)₂`.
*   To find the structure of Compound A, we perform a "retro-Wittig" analysis: we replace the `=C(CH₃)₂` group in the product with a carbonyl group (`=O`).
*   This gives us a ketone at position 4 of a 1,2-dimethylcyclopentane ring. According to IUPAC nomenclature, the carbonyl group is given the highest priority and numbered as position 1. This places the methyl groups at positions 3 and 4.
*   Therefore, **Compound A is 3,4-dimethylcyclopentan-1-one**.
*   **Hint (b)** confirms this. The IR peak for Compound A is at ~1750 cm⁻¹, which is characteristic of a five-membered ring ketone (a cyclopentanone) due to increased ring strain.

**2. Tracing the Reactions (A → E):**
*   **A → B:** Compound A (3,4-dimethylcyclopentan-1-one) reacts with HCN to form a cyanohydrin. **Compound B is 1-cyano-3,4-dimethylcyclopentan-1-ol**.
*   **B → C:** The nitrile group (-CN) of Compound B is reduced with H₂/Pd to a primary amine (-CH₂NH₂). **Compound C is 1-(aminomethyl)-3,4-dimethylcyclopentan-1-ol**.
*   **C → D → E:** This is the Tiffeneau-Demjanov rearrangement.
    *   The primary amine (Compound C) reacts with nitrous acid (HNO₂) to form an unstable diazonium salt (Compound D).
    *   Compound D loses nitrogen gas (N₂, the "inert gas") to form a highly unstable primary carbocation.
    *   This carbocation immediately rearranges. A carbon-carbon bond from the cyclopentane ring migrates, expanding the five-membered ring into a more stable six-membered ring.
    *   The final product is a ketone where the methyl groups, originally at positions 3 and 4, are now at positions 3 and 4 of the new six-membered ring.
*   Therefore, **Compound E is 3,4-dimethylcyclohexan-1-one**.

**3. Final Verification:**
*   **Hint (b)** states that the IR peak for Compound E is at ~1715 cm⁻¹. This frequency is characteristic of a six-membered ring ketone (a cyclohexanone), which has less ring strain than a cyclopentanone. This spectroscopic data strongly confirms that a ring expansion from a five-membered ring (A) to a six-membered ring (E) has occurred.

### **Conclusion from Candidate Answers**

A review of the provided candidate answers shows an overwhelming consensus on the chemical pathway. Nearly all analyses correctly identify Compound A as 3,4-dimethylcyclopentan-1-one and Compound E as 3,4-dimethylcyclohexan-1-one through the Tiffeneau-Demjanov rearrangement.

While several candidate answers have conclusions (e.g., <<<D>>>, <<<A>>>, <<<C>>>) that contradict their own detailed reasoning, this appears to be due to errors in mapping the correct chemical name to the final letter choice. The chemical logic presented consistently points to one structure.

The correct structure, 3,4-dimethylcyclohexan-1-one, corresponds to option B.

<<<B>>>
"""

# Instantiate and run the checker
checker = ChemistryProblemChecker(problem_details, llm_final_answer)
result = checker.check()
print(result)