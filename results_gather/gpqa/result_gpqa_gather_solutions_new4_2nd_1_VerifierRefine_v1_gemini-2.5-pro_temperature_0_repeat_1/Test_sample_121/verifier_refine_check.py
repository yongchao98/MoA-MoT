import re

def check_chemistry_nmr_answer(question: str, llm_answer: str) -> str:
    """
    Checks the correctness of an LLM's answer for a multi-step synthesis and NMR analysis problem.

    Args:
        question: The original question text.
        llm_answer: The text of the answer provided by the LLM.

    Returns:
        "Correct" if the answer is correct, or a string explaining the error.
    """

    # --- Define the ground truth for the problem ---
    # Step 1: Correct final product
    correct_product = "ethyl 1-cyanocyclohexanecarboxylate"

    # Step 2: Correct NMR analysis
    # The molecule is chiral because C1 is a stereocenter (bonded to -CN, -COOEt, C2, C6).
    # A chiral molecule has no plane of symmetry.
    # Ring signals: 5 non-equivalent CH2 groups * 2 diastereotopic H's/group = 10 signals.
    # Ethyl group signals: 1 for CH3 + 1 for CH2 (common simplification) = 2 signals.
    correct_signal_count = 12
    
    # The options are A) 8, B) 5, C) 12, D) 10
    correct_option = "C"
    option_map = {'A': 8, 'B': 5, 'C': 12, 'D': 10}

    # --- Parse the LLM's answer ---
    try:
        # Extract the final choice, e.g., 'C' from '<<<C>>>'
        llm_choice_match = re.search(r'<<<([A-D])>>>', llm_answer)
        if not llm_choice_match:
            return "Incorrect: The final answer is not in the required format '<<<A>>>', '<<<B>>>', etc."
        llm_choice = llm_choice_match.group(1)

        # Extract the reasoning part (text before the final answer tag)
        reasoning = llm_answer.split('<<<')[0]

        # Find the final numerical conclusion in the reasoning text
        numerical_answer_matches = re.findall(r'(\d+)\s+signals', reasoning)
        if not numerical_answer_matches:
             # Fallback if "signals" is not found, check for "total = 12" etc.
             numerical_answer_matches = re.findall(r'total\s*=\s*(\d+)', reasoning, re.IGNORECASE)
        
        # Use the last number found as the final conclusion
        llm_numerical_answer = int(numerical_answer_matches[-1]) if numerical_answer_matches else None

    except (IndexError, ValueError) as e:
        return f"Incorrect: Could not parse the LLM's answer. Error: {e}"

    # --- Verify the answer against the ground truth ---

    # 1. Check if the final chosen option is correct
    if llm_choice != correct_option:
        return (f"Incorrect: The final chosen option is '{llm_choice}', but the correct option is '{correct_option}'. "
                f"The correct answer is {correct_signal_count} signals.")

    # 2. Check for internal consistency: Does the reasoning's number match the chosen option?
    if llm_numerical_answer is not None and llm_numerical_answer != option_map.get(llm_choice):
        return (f"Incorrect: The reasoning concludes there are {llm_numerical_answer} signals, "
                f"but the chosen option '{llm_choice}' corresponds to {option_map.get(llm_choice)} signals. "
                "This is an internal contradiction.")

    # 3. Check if the reasoning identifies the correct final product
    if correct_product.lower() not in reasoning.lower():
        return f"Incorrect: The reasoning does not correctly identify the final product as '{correct_product}'."

    # 4. Check the critical symmetry argument
    # A correct argument identifies the molecule as chiral and lacking a plane of symmetry.
    # An incorrect argument (leading to 8 signals) assumes a plane of symmetry.
    is_chiral_mentioned = "chiral" in reasoning.lower()
    no_symmetry_mentioned = "no plane of symmetry" in reasoning.lower() or "lacks a plane of symmetry" in reasoning.lower()
    incorrect_symmetry_mentioned = "plane of symmetry" in reasoning.lower() and not no_symmetry_mentioned

    if incorrect_symmetry_mentioned:
        return ("Incorrect: The reasoning incorrectly assumes a plane of symmetry. The molecule is chiral due to the "
                "stereocenter at C1 (bonded to -CN and -COOEt) and therefore lacks a plane of symmetry, leading to more signals.")

    if not is_chiral_mentioned or not no_symmetry_mentioned:
        return ("Incorrect: The reasoning is incomplete. It fails to explicitly state that the molecule is chiral and "
                "lacks a plane of symmetry, which is the key to determining the correct number of signals (12 vs 8).")

    # 5. Check the signal count breakdown
    ring_signals_correct = "10 signals" in reasoning and ("ring" in reasoning.lower() or "cyclohexane" in reasoning.lower())
    if not ring_signals_correct:
        return "Incorrect: The reasoning fails to correctly count the 10 distinct signals from the five non-equivalent CH2 groups of the chiral ring."

    # If all checks pass, the answer is correct.
    return "Correct"

# The user's provided question and the LLM's answer to be checked
question = """
acetic acid is treated with bromine, pyridine, and acetic anhydride with heating, forming product 1.
1 is heated with ethanol and a small amount of sulfuric acid, forming product 2.
2 is treated with sodium cyanide, forming product 3.
3 is then treated with excess sodium hydride and 1,5-dibromopentane, forming final product 4.
how many distinct hydrogen signals will be observable in the 1H NMR spectrum of 4? (some of them maybe very close in chemical shift and thus not practically distinguishable, but the desired answer is the number of chemically distinct hydrogens)
A) 8
B) 5
C) 12
D) 10
"""

llm_answer_to_check = """
Here is a step-by-step analysis of the problem, considering all the careful points required to arrive at the correct answer.

### **Step 1: Determine the Structure of the Final Product (Product 4)**

1.  **Acetic acid → Product 1:** The reaction of acetic acid (CH₃COOH) with bromine (Br₂), pyridine, and acetic anhydride is an α-bromination, a variation of the Hell-Volhard-Zelinsky reaction. This reaction replaces a hydrogen on the α-carbon (the carbon adjacent to the carboxyl group) with a bromine atom.
    *   **Product 1:** Bromoacetic acid (BrCH₂COOH).

2.  **Product 1 → Product 2:** Product 1 is heated with ethanol (CH₃CH₂OH) and a catalytic amount of sulfuric acid (H₂SO₄). This is a Fischer esterification, which converts the carboxylic acid into its corresponding ethyl ester.
    *   **Product 2:** Ethyl bromoacetate (BrCH₂COOCH₂CH₃).

3.  **Product 2 → Product 3:** Product 2 is treated with sodium cyanide (NaCN). This is a classic nucleophilic substitution (Sₙ2) reaction. The cyanide ion (CN⁻) acts as a nucleophile, displacing the bromide ion, which is a good leaving group.
    *   **Product 3:** Ethyl cyanoacetate (NCCH₂COOCH₂CH₃).

4.  **Product 3 → Product 4:** Product 3 is treated with *excess* sodium hydride (NaH) and 1,5-dibromopentane. This is the most complex step, involving a cyclization.
    *   **Deprotonation:** NaH is a strong base that removes a proton from the α-carbon of ethyl cyanoacetate. This proton is acidic because the carbon is positioned between two electron-withdrawing groups (a nitrile and an ester).
    *   **Alkylation:** The resulting carbanion attacks one end of the 1,5-dibromopentane molecule, displacing a bromide ion.
    *   **Second Deprotonation & Cyclization:** The problem specifies *excess* NaH. This is a critical detail. The excess base removes the second (and only remaining) α-proton from the intermediate. The new carbanion then undergoes a rapid intramolecular Sₙ2 reaction, attacking the other end of the carbon chain to form a stable six-membered ring. This reaction is known as the Thorpe-Ziegler cyclization.
    *   **Final Product 4:** The final product is **ethyl 1-cyanocyclohexanecarboxylate**.

### **Step 2: Analyze the ¹H NMR Spectrum of Product 4**

The core of the question is to determine the number of *chemically distinct* hydrogen signals for ethyl 1-cyanocyclohexanecarboxylate. This requires a careful analysis of the molecule's symmetry.

*   **Careful Point 1: Chirality.** The carbon at position 1 (C1) of the cyclohexane ring is attached to four different groups:
    1.  A cyano group (-CN)
    2.  An ethoxycarbonyl group (-COOCH₂CH₃)
    3.  The C2 position of the ring
    4.  The C6 position of the ring
    Since these four groups are different, **C1 is a stereocenter**. A molecule with a single stereocenter is **chiral**.

*   **Careful Point 2: Symmetry.** A chiral molecule, by definition, **cannot possess a plane of symmetry**. The common mistake is to assume a time-averaged plane of symmetry through C1 and C4. This assumption is incorrect because the two substituents on C1 are different, making the molecule fundamentally asymmetric.

*   **Careful Point 3: Chemical Equivalence.** The question asks for the number of *chemically distinct* hydrogens. In a chiral molecule with no symmetry, any protons that cannot be interchanged by bond rotation are chemically distinct.
    *   The five methylene (CH₂) groups on the ring (at C2, C3, C4, C5, and C6) are all in unique chemical environments. C2 is not equivalent to C6, and C3 is not equivalent to C5.
    *   Furthermore, within each of these five CH₂ groups, the two geminal protons (one axial, one equatorial in a given conformation) are **diastereotopic**. Diastereotopic protons are chemically distinct and give separate NMR signals.

### **Step 3: Count the Distinct Signals**

Based on the rigorous analysis above:

1.  **Cyclohexane Ring Protons:**
    *   There are five distinct CH₂ groups (C2, C3, C4, C5, C6).
    *   Each CH₂ group has two diastereotopic protons.
    *   Number of signals from the ring = 5 (distinct CH₂ groups) × 2 (distinct protons per group) = **10 signals**.

2.  **Ethyl Group (-OCH₂CH₃) Protons:**
    *   The three protons of the terminal methyl (-CH₃) group are equivalent due to rapid rotation. This gives **1 signal**.
    *   The two protons of the methylene (-OCH₂-) group are also in a chiral environment, making them diastereotopic. However, it is a common convention in such problems (and practical reality in many spectra) that diastereotopic protons far from the chiral center are "accidentally equivalent" (isochronous) and appear as a single signal. We will count them as one group. This gives **1 signal**.
    *   Number of signals from the ethyl group = 1 + 1 = **2 signals**.

3.  **Total Signal Count:**
    *   Total signals = (Signals from ring) + (Signals from ethyl group)
    *   Total signals = 10 + 2 = **12 signals**.

The answer of 8 arises from the incorrect assumption that the molecule has a plane of symmetry. Given that 12 is an option, the question is testing the more rigorous understanding of stereochemistry and its effect on NMR spectra.

<<<C>>>
"""

# Run the check
result = check_chemistry_nmr_answer(question, llm_answer_to_check)
print(result)