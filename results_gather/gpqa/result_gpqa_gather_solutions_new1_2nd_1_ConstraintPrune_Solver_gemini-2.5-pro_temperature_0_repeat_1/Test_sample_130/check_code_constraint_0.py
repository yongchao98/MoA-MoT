import re

def check_correctness_of_chemistry_answer(question: str, llm_answer: str) -> str:
    """
    Checks the correctness of an LLM's answer to a complex organic chemistry problem.

    The function codifies the key chemical principles needed to solve the problem:
    1.  Reactant identification.
    2.  Stereoselectivity of the Diels-Alder reaction, considering steric hindrance.
    3.  Assignment of NMR signals.
    4.  Prediction of NOESY correlations based on the 3D structure of the major product.

    Args:
        question: The original question text.
        llm_answer: The text of the LLM's answer, including the final choice.

    Returns:
        A string indicating "Correct" or a detailed reason for the incorrectness.
    """

    # --- Step 1: Define the chemical knowledge base for this specific problem ---

    # Proton signal assignments based on typical chemical shifts for the product structure
    signal_assignments = {
        "H_anhydride": "2H singlet at ~3.5 ppm",
        "Me_vinyl": "6H singlet at ~1.7 ppm",
        "H_bridge": "1H doublet at ~1.5 ppm",
        "Me_bridgehead": "6H singlet at ~1.0 ppm"
    }

    # Define the answer options by mapping the letter to the set of interacting signals
    # This mapping is based on the options provided in the user's prompt.
    options = {
        "A": {signal_assignments["Me_vinyl"], signal_assignments["H_anhydride"]},
        "B": {signal_assignments["Me_bridgehead"], signal_assignments["H_bridge"]},
        "C": {signal_assignments["H_bridge"], signal_assignments["H_anhydride"]},
        "D": {signal_assignments["Me_bridgehead"], signal_assignments["Me_vinyl"]}
    }

    # --- Step 2: Extract the final answer from the LLM's response ---
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "Failure: Could not find a final answer in the format <<<A>>>, <<<B>>>, etc. in the provided text."
    
    final_answer_letter = match.group(1)

    # --- Step 3: Apply chemical principles to determine the correct answer ---

    # Principle 1: The diene (1,2,3,4-tetramethyl-1,3-cyclopentadiene) is extremely bulky.
    # Principle 2: Due to severe steric hindrance, the Alder-endo rule is overridden.
    # The major product of this specific Diels-Alder reaction is the EXO adduct.
    major_product_isomer = "exo"

    # Principle 3: Determine the key NOESY interaction for the major (exo) product.
    # In the EXO adduct, the anhydride protons (H_anhydride) are on the endo face,
    # spatially close to the vinylic methyl groups (Me_vinyl), which are also on the endo face.
    # The characteristic interaction for the ENDO adduct would be between H_anhydride and H_bridge.
    
    if major_product_isomer == "exo":
        # Protons interacting in the major product
        interacting_protons = {"Me_vinyl", "H_anhydride"}
    else: # This case is not reached based on correct chemistry, but included for completeness
        interacting_protons = {"H_bridge", "H_anhydride"}

    # Translate the interacting protons to their corresponding signal descriptions
    correct_signals = {signal_assignments[proton] for proton in interacting_protons}

    # --- Step 4: Compare the LLM's answer with the derived correct answer ---
    
    llm_chosen_signals = options.get(final_answer_letter)

    if llm_chosen_signals == correct_signals:
        return "Correct"
    else:
        # Find which option letter *should* have been correct
        correct_letter = "Unknown"
        for letter, signal_set in options.items():
            if signal_set == correct_signals:
                correct_letter = letter
                break
        
        reason = (
            f"Incorrect. The provided answer is '{final_answer_letter}', but the correct answer is '{correct_letter}'.\n"
            f"Reasoning:\n"
            f"1. The most critical point is determining the major product. While the Alder-endo rule is common, it is overridden here by the severe steric hindrance of the 1,2,3,4-tetramethyl-1,3-cyclopentadiene diene. Therefore, the **exo adduct is the major product**.\n"
            f"2. In the 3D structure of the major (exo) product, the anhydride protons (~3.5 ppm) and the vinylic methyl protons (~1.7 ppm) are on the same face of the molecule and are spatially close.\n"
            f"3. This proximity results in a NOESY cross-peak between the '6H singlet at ~1.7 ppm' and the '2H singlet at ~3.5 ppm', which corresponds to option '{correct_letter}'.\n"
            f"4. The chosen answer '{final_answer_letter}' describes an interaction between {', '.join(llm_chosen_signals)}, which is incorrect for the major product."
        )
        return reason

# The user's provided answer to check
final_answer_from_prompt = """
Here is a step-by-step analysis of the problem, synthesizing the correct reasoning from the provided answers and identifying critical points.

### Step-by-step Analysis

1.  **Identify the Reactants:**
    *   **Dienophile:** The problem describes the dehydration of a *cis*-alkene dicarboxylic acid to an anhydride with specific NMR data (¹H: single peak at 7 ppm; ¹³C: two peaks at 137 ppm and 165 ppm). The high symmetry and characteristic chemical shifts unambiguously identify this compound as **maleic anhydride**.
    *   **Diene:** The diene is explicitly named as **1,2,3,4-tetramethyl-1,3-cyclopentadiene**.

2.  **Identify the Reaction and Determine the Major Product:**
    *   The reaction is a **Diels-Alder [4+2] cycloaddition**.
    *   This reaction produces two stereoisomers: an *endo* and an *exo* adduct. A crucial step is to determine which is the major product.
    *   While the **Alder-endo rule** typically predicts the *endo* adduct as the major product due to favorable electronic interactions, this rule can be overridden by severe steric hindrance.
    *   The diene, 1,2,3,4-tetramethyl-1,3-cyclopentadiene, is exceptionally bulky. The four methyl groups create significant steric clash in the *endo* transition state. The *exo* transition state, where the anhydride approaches from the less hindered face, is much more favorable.
    *   Therefore, due to this overwhelming steric hindrance, the **exo adduct is the major product**, and the *endo* adduct is the minor product.

3.  **Analyze the 3D Structure and Predict the NOESY Correlation:**
    *   A NOESY cross-peak indicates that two protons are close in space (< 5 Å). We are looking for a cross-peak that is present in the major (*exo*) product but absent or weak in the minor (*endo*) product.
    *   **Signal Assignments:**
        *   **~3.5 ppm (2H singlet):** The two equivalent protons on the anhydride ring (`H_anhydride`).
        *   **~1.7 ppm (6H singlet):** The six protons of the two equivalent methyl groups on the double bond (`Me_vinyl`).
        *   **~1.5 ppm (1H doublet):** One of the two non-equivalent protons on the C7 bridge (`H_bridge`).
        *   **~1.0 ppm (6H singlet):** The six protons of the two equivalent methyl groups at the bridgehead positions (`Me_bridgehead`).
    *   **Spatial Analysis of the Major (*Exo*) Product:** In the *exo* adduct, the anhydride ring is positioned *anti* to the C7 bridge. This geometry places the `H_anhydride` protons on the *endo* face of the bicyclic system, pointing towards the double bond. The `Me_vinyl` groups are also on the *endo* face. This brings them into close spatial proximity.
    *   **Spatial Analysis of the Minor (*Endo*) Product:** In the *endo* adduct, the `H_anhydride` protons are on the *exo* face, far away from the `Me_vinyl` groups.

4.  **Conclusion:**
    *   The unique NOESY cross-peak in the major product's spectrum must arise from the interaction between the anhydride protons and the vinylic methyl protons.
    *   This corresponds to a cross-peak connecting the **2H singlet at ~3.5 ppm** and the **6H singlet at ~1.7 ppm**.
    *   This interaction is described in option A.

### Careful Points in the Analysis

1.  **Reactant Identification:** The first step is to correctly identify the reactants. The NMR data for the anhydride strongly points to maleic anhydride, which is a common and correct starting point for all candidate answers.
2.  **Stereoselectivity (The Most Critical Point):** The central challenge of this problem is determining the major product. A common mistake is to automatically apply the Alder-endo rule. It is essential to recognize that the extreme steric bulk of the 1,2,3,4-tetramethyl-1,3-cyclopentadiene diene overrides the endo rule, making the **exo adduct the major product**. Any analysis based on the *endo* adduct being major is fundamentally flawed.
3.  **3D Visualization and NOE Prediction:** Correctly visualizing the 3D structure of the major (*exo*) adduct is necessary to predict the NOESY correlation. The key insight is that in this isomer, the anhydride protons (`H_anhydride`) and the vinylic methyl protons (`Me_vinyl`) are on the same face of the molecule and therefore spatially close.
4.  **Matching Interaction to Options:** The final step is to accurately match the identified interacting protons and their corresponding NMR signals to the given options. The interaction is between the 6H singlet at ~1.7 ppm (`Me_vinyl`) and the 2H singlet at ~3.5 ppm (`H_anhydride`), which corresponds to option A. Several candidate answers correctly identified the interaction but failed to match it to the correct lettered option.

<<<A>>>
"""

# Run the check
result = check_correctness_of_chemistry_answer("The chemistry question", final_answer_from_prompt)
print(result)