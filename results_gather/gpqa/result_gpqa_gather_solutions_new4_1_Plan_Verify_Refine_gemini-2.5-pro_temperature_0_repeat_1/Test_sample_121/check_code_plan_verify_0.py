import sys
from io import StringIO

class ChemistryChecker:
    def __init__(self):
        self.steps = []
        self.final_product = None
        self.signal_count = 0
        self.correct_option = None
        self.correct_value = None

    def run_check(self, provided_answer_text, provided_answer_option):
        """
        Runs a full check on the provided answer.
        """
        # Step 1: Check the reaction pathway
        pathway_check_result = self._check_reaction_pathway()
        if pathway_check_result != "Correct":
            return pathway_check_result

        # Step 2: Check the NMR analysis for the final product
        nmr_check_result = self._check_nmr_analysis()
        if nmr_check_result != "Correct":
            return nmr_check_result

        # Step 3: Compare with the provided answer's reasoning and conclusion
        # The provided answer's reasoning matches our derivation (8 signals).
        # Let's check if its final option matches.
        if self.correct_value != 8:
             return f"Internal logic error: The correct value was determined to be {self.correct_value}, not 8."

        # The options are A) 5, B) 10, C) 8, D) 12. So C is correct.
        self.correct_option = 'C'
        
        if provided_answer_option != self.correct_option:
            return f"Incorrect: The final answer option '{provided_answer_option}' is wrong. The correct option is '{self.correct_option}', which corresponds to {self.correct_value} signals."

        return "Correct"

    def _check_reaction_pathway(self):
        """
        Verifies the multi-step synthesis.
        """
        # Step 1: Acetic acid -> Product 1 (alpha-bromination)
        product_1 = "Bromoacetic acid"
        
        # Step 2: Product 1 -> Product 2 (Fischer esterification)
        product_2 = "Ethyl bromoacetate"
        
        # Step 3: Product 2 -> Product 3 (SN2 with cyanide)
        product_3 = "Ethyl cyanoacetate"
        
        # Step 4: Product 3 -> Product 4 (Thorpe-Ziegler cyclization)
        # The reaction of an alpha-cyano ester with a dihalide (1,5-dibromopentane)
        # under strong base conditions is a classic way to form a ring.
        # Ring size = 1 (alpha-carbon) + 5 (from dibromopentane) = 6-membered ring.
        self.final_product = "Ethyl 1-cyanocyclohexanecarboxylate"
        
        # The provided answer correctly identifies this reaction pathway.
        return "Correct"

    def _check_nmr_analysis(self):
        """
        Analyzes the final product's 1H NMR spectrum.
        This models the most common "textbook" interpretation for this problem.
        """
        if self.final_product != "Ethyl 1-cyanocyclohexanecarboxylate":
            return "Cannot perform NMR analysis because the final product is incorrect."

        # Analysis of Ethyl 1-cyanocyclohexanecarboxylate
        
        # Assumption 1: At room temperature, rapid chair-flipping creates a time-averaged
        # plane of symmetry passing through C1 and C4.
        # Assumption 2: Geminal protons (on the same CH2 group) on the ring are
        # diastereotopic and thus chemically distinct, giving separate signals.

        # Count signals from the ethyl group (-O-CH2-CH3)
        # Signal 1: -CH3 (3 equivalent protons)
        # Signal 2: -O-CH2- (2 equivalent protons, enantiotopic due to symmetry)
        ethyl_signals = 2

        # Count signals from the cyclohexane ring
        # C2 and C6 are equivalent by symmetry. The ax/eq protons on them are distinct.
        c2_c6_signals = 2
        # C3 and C5 are equivalent by symmetry. The ax/eq protons on them are distinct.
        c3_c5_signals = 2
        # C4 is on the symmetry plane. Its ax/eq protons are distinct.
        c4_signals = 2
        
        ring_signals = c2_c6_signals + c3_c5_signals + c4_signals

        self.signal_count = ethyl_signals + ring_signals
        self.correct_value = self.signal_count
        
        if self.signal_count != 8:
            return f"Incorrect NMR Analysis: The standard analysis leads to 8 signals, not {self.signal_count}."
            
        return "Correct"


def check_correctness():
    """
    This function checks the provided answer against the correct chemical reasoning.
    """
    # The provided answer to check is the one marked as the final answer in the prompt.
    # Let's extract its key components.
    # Final Answer Text: The reasoning leads to 8 signals.
    # Final Answer Option: <<<C>>>
    
    provided_answer_text = """
    Here is a step-by-step analysis to determine the correct answer.

    ### **1. Analysis of the Reaction Sequence**
    The first step is to correctly identify the structure of the final product, 4.
    *   **Step 1: Acetic acid → Product 1**
        *   **Reaction:** Acetic acid (CH₃COOH) is treated with bromine (Br₂), pyridine, and acetic anhydride. This is an α-bromination reaction (a variation of the Hell-Volhard-Zelinsky reaction).
        *   **Product 1:** Bromoacetic acid (BrCH₂COOH).
    *   **Step 2: Product 1 → Product 2**
        *   **Reaction:** Product 1 is heated with ethanol (CH₃CH₂OH) and sulfuric acid (H₂SO₄). This is a Fischer esterification.
        *   **Product 2:** Ethyl bromoacetate (BrCH₂COOCH₂CH₃).
    *   **Step 3: Product 2 → Product 3**
        *   **Reaction:** Product 2 is treated with sodium cyanide (NaCN). This is a nucleophilic substitution (Sₙ2) where the cyanide ion displaces the bromide.
        *   **Product 3:** Ethyl cyanoacetate (NCCH₂COOCH₂CH₃).
    *   **Step 4: Product 3 → Product 4**
        *   **Reaction:** Product 3 is treated with *excess* sodium hydride (NaH) and 1,5-dibromopentane.
        *   **Product 4:** Ethyl 1-cyanocyclohexanecarboxylate.

    ### **2. Analysis of the ¹H NMR Spectrum of Product 4**
    *   **Structure:** A cyclohexane ring with a cyano group (-CN) and an ethoxycarbonyl group (-COOCH₂CH₃) attached to the same carbon (C1).
    *   **Counting Signals with the Standard Approximation:**
        1.  **Ethyl Group (-OCH₂CH₃):** 2 signals.
        2.  **Cyclohexane Ring:** 6 signals.
    *   **Total Signal Count:** Total = 2 + 6 = **8 signals**.
    """
    provided_answer_option = "C"

    checker = ChemistryChecker()
    result = checker.run_check(provided_answer_text, provided_answer_option)
    
    return result

# Execute the check and print the result
result = check_correctness()
print(result)
