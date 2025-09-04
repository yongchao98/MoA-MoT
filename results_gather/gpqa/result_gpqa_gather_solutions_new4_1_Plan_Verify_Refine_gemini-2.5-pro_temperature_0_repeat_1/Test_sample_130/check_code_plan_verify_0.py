import re

class DielsAlderNOESYChecker:
    """
    A class to verify the answer to a complex stereochemistry and NMR problem.
    It encodes the chemical principles needed to solve the problem and compares
    the expected outcome with the provided answer.
    """

    def __init__(self, question, llm_answer_text):
        self.question = question
        self.llm_answer_text = llm_answer_text
        self.final_choice = self._parse_final_choice(llm_answer_text)

        # Step 1: Define proton signals and map them to chemical groups.
        self.signal_assignments = {
            'H_anhyd': '2H singlet at ~3.5 ppm',
            'Me_vinyl': '6H singlet at ~1.7 ppm',
            'Me_bridgehead': '6H singlet at ~1.0 ppm',
            'H_bridge': '1H doublet at ~1.5 ppm'
        }

        # Map the options from the question to the chemical groups involved.
        self.options_map = {
            'A': {'Me_bridgehead', 'Me_vinyl'},
            'B': {'H_bridge', 'H_anhyd'},
            'C': {'Me_bridgehead', 'H_bridge'},
            'D': {'Me_vinyl', 'H_anhyd'}
        }

        # Step 2: Encode the spatial proximity rules for NOESY in each isomer.
        self.noe_correlations = {
            'endo': {'protons': {'H_anhyd', 'H_bridge'},
                     'reason': 'In the endo isomer, the anhydride protons are spatially close to the C7 bridge protons.'},
            'exo': {'protons': {'H_anhyd', 'Me_vinyl'},
                    'reason': 'In the exo isomer, the anhydride protons are spatially close to the vinylic methyl groups.'}
        }

    def _parse_final_choice(self, text):
        """Extracts the final answer choice (e.g., 'D') from the text."""
        match = re.search(r'<<<([A-D])>>>', text)
        return match.group(1) if match else None

    def check_correctness(self):
        """
        Runs the verification logic to check if the provided answer is correct.
        """
        if not self.final_choice:
            return "Incorrect. The final answer is not provided in the required format '<<<X>>>'."

        # Step 3: Determine the major product based on chemical principles.
        # The standard "endo rule" is reversed due to severe steric hindrance from the
        # four methyl groups on the diene.
        major_product_type = 'exo'
        
        # Determine the expected NOESY cross-peak for the major product.
        expected_noe = self.noe_correlations[major_product_type]
        expected_protons = expected_noe['protons']

        # Get the protons corresponding to the answer's chosen option.
        chosen_protons = self.options_map.get(self.final_choice)

        # Step 4: Compare the expected result with the answer's choice.
        if chosen_protons == expected_protons:
            return "Correct"
        else:
            # Construct a detailed reason for the incorrectness.
            correct_option = None
            for option, protons in self.options_map.items():
                if protons == expected_protons:
                    correct_option = option
                    break
            
            reason = (
                f"The answer is incorrect.\n"
                f"1. Major Product Determination: Due to severe steric hindrance from the tetramethyl-substituted diene, the major product is the '{major_product_type}' adduct, not the endo adduct.\n"
                f"2. Expected NOESY Correlation: {expected_noe['reason']}\n"
                f"3. Conclusion: The unique cross-peak in the major product should connect the anhydride protons (~3.5 ppm) and the vinylic methyl protons (~1.7 ppm). This corresponds to option {correct_option}.\n"
                f"The provided answer chose option {self.final_choice}, which is incorrect."
            )
            return reason

# --- Execution ---
# The user's provided question and the LLM's answer to be checked.
question_text = """
A dicarboxylic acid containing a cis-alkene was dehydrated to the corresponding anhydride. This anhydride, which has a single peak (7 ppm) in its 1H NMR spectrum and two peaks (137 ppm, 165 ppm) in its 13C spectrum, was reacted with 1,2,3,4-tetramethyl-1,3-cyclopentadiene to yield two products (a major and a minor product). Both products have a molecular formula of C13H16O3.

The 2D NOESY NMR spectra of the major and minor products were taken. A cross-peak was present in the spectrum of the major product but absent (or very weak) in spectrum of the minor product.

In the spectrum of the major product, which two resonances are connected by this cross-peak (integrals are taken from the corresponding 1H NMR)?

A) A 6H singlet at ~1 ppm and a 6H singlet at ~1.7 ppm
B) A 1H doublet at ~1.5 ppm and a 2H singlet at ~3.5 ppm
C) A 6H singlet at ~1 ppm and a 1H doublet at ~1.5 ppm
D) A 6H singlet at ~1.7 ppm and a 2H singlet at ~3.5 ppm
"""

llm_answer_to_check = """
Here is a step-by-step analysis of the problem to determine the final answer.

### Step 1: Identify the Reactants
*   **Dienophile:** The problem describes a *cis*-dicarboxylic acid that is dehydrated to an anhydride. The NMR data for this anhydride (¹H: single peak at 7 ppm; ¹³C: two peaks at 137 ppm and 165 ppm) are characteristic of **maleic anhydride**. The high symmetry of maleic anhydride explains why there is only one type of proton and two types of carbons.
*   **Diene:** The diene is explicitly named as **1,2,3,4-tetramethyl-1,3-cyclopentadiene**.

### Step 2: The Diels-Alder Reaction and Stereoselectivity
*   The reaction is a [4+2] cycloaddition between the diene and dienophile, forming a bicyclo[2.2.1]heptene derivative.
*   Two stereoisomers are formed: the *endo* and *exo* adducts. The question is which one is the major product.
*   The standard "Alder-endo rule" predicts the *endo* isomer to be the major kinetic product due to favorable secondary orbital interactions.
*   However, this rule can be overridden by severe steric hindrance. The diene, 1,2,3,4-tetramethyl-1,3-cyclopentadiene, is extremely bulky. An *endo* approach would force the anhydride ring to pass by the methyl groups on the diene, creating significant steric clash. The *exo* approach is much less sterically hindered.
*   Therefore, in this specific case, the **exo adduct is the major product**, and the *endo* adduct is the minor product.

### Step 3: Analyze the 3D Structures and NOESY Correlations
A NOESY cross-peak indicates that two protons are close in space (< 5 Å). We are looking for an interaction present in the major (*exo*) product but absent in the minor (*endo*) product.

*   **Proton Signal Assignments:** Based on typical chemical shifts for this type of structure:
    *   **2H singlet at ~3.5 ppm:** The two equivalent protons on the anhydride ring (H_anhyd).
    *   **6H singlet at ~1.7 ppm:** The six protons of the two equivalent methyl groups on the double bond (Me_vinyl).
    *   **6H singlet at ~1.0 ppm:** The six protons of the two equivalent methyl groups at the bridgehead positions (Me_bridgehead).
    *   **1H doublet at ~1.5 ppm:** One of the two non-equivalent protons on the C7 bridge (H_bridge).

*   **Spatial Proximity in the *Exo* Isomer (Major Product):**
    *   In the *exo* structure, the anhydride ring is positioned on the same side as the C7 bridge, pointing away from the C=C double bond.
    *   This geometry brings the anhydride protons (**H_anhyd**, ~3.5 ppm) very close to the vinylic methyl groups (**Me_vinyl**, ~1.7 ppm).
    *   A strong NOESY cross-peak is expected between these two signals.

*   **Spatial Proximity in the *Endo* Isomer (Minor Product):**
    *   In the *endo* structure, the anhydride ring is tucked underneath the C=C double bond.
    *   This places the anhydride protons (**H_anhyd**, ~3.5 ppm) far from the vinylic methyl groups.
    *   Instead, they are close to one of the C7 bridge protons (**H_bridge**, ~1.5 ppm). A cross-peak would be expected here for the minor product.

### Step 4: Conclusion
The question asks for the cross-peak present in the spectrum of the **major product**. Our analysis shows:
1.  The major product is the **exo** adduct.
2.  The characteristic NOESY cross-peak for the *exo* adduct is between the anhydride protons and the vinylic methyl protons.
3.  These protons correspond to the **2H singlet at ~3.5 ppm** and the **6H singlet at ~1.7 ppm**.

This matches option D.

<<<D>>>
"""

# Instantiate the checker and run the verification.
checker = DielsAlderNOESYChecker(question_text, llm_answer_to_check)
result = checker.check_correctness()
print(result)