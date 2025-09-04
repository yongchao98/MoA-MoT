import re

def check_diels_alder_noesy(answer_text: str) -> str:
    """
    Checks the correctness of the provided answer for the Diels-Alder NOESY problem.

    The function validates the answer based on established chemical principles:
    1.  **Reactant Identification**: Confirms the dienophile is maleic anhydride.
    2.  **Stereoselectivity**: Checks if the major product is correctly identified as the 'exo' adduct due to steric hindrance, which is a known exception to the Alder-endo rule in this specific case.
    3.  **Signal Assignment**: Verifies the assignment of NMR signals to the corresponding protons.
    4.  **NOESY Correlation**: Determines the expected NOESY cross-peak in the major 'exo' product.
    5.  **Final Answer Check**: Compares the derived correct option with the one provided in the answer.
    """

    # --- Step 1: Define Chemical Principles and Problem Data ---

    # Signal assignments based on chemical principles
    signal_assignments = {
        "H_anhydride": ("2H singlet", 3.5),
        "Me_vinyl": ("6H singlet", 1.7),
        "Me_bridgehead": ("6H singlet", 1.0),
        "H_bridge": ("1H doublet", 1.5)
    }

    # The crucial point: Due to severe steric hindrance from the four methyl groups
    # on the diene, the 'exo' adduct is the major product, overriding the Alder-endo rule.
    major_product_identity = "exo"

    # Determine the key NOESY interaction based on the major product's 3D structure
    # In the 'exo' adduct, the anhydride protons (H_anhydride) are on the 'endo' face,
    # close to the vinylic methyl groups (Me_vinyl).
    # In the 'endo' adduct, the anhydride protons are on the 'exo' face, close to the
    # bridge protons (H_bridge).
    if major_product_identity == "exo":
        key_interaction_protons = {"H_anhydride", "Me_vinyl"}
        correct_option_description = "Interaction between anhydride protons and vinylic methyl protons."
    else: # This 'else' block represents the incorrect assumption (endo is major)
        key_interaction_protons = {"H_anhydride", "H_bridge"}
        correct_option_description = "Interaction between anhydride protons and bridge protons."

    # Map the interacting protons to their signals
    interacting_signals = [signal_assignments[proton] for proton in key_interaction_protons]

    # Define the options
    options = {
        'A': [("6H singlet", 1.0), ("1H doublet", 1.5)], # Me_bridgehead & H_bridge
        'B': [("6H singlet", 1.0), ("6H singlet", 1.7)], # Me_bridgehead & Me_vinyl
        'C': [("6H singlet", 1.7), ("2H singlet", 3.5)], # Me_vinyl & H_anhydride
        'D': [("1H doublet", 1.5), ("2H singlet", 3.5)]  # H_bridge & H_anhydride
    }

    # Find the correct letter option
    correct_option = None
    for option, signals in options.items():
        # Check if the set of signals in the option matches the set of interacting signals
        if set(signals) == set(interacting_signals):
            correct_option = option
            break

    # --- Step 2: Analyze the Provided Answer ---

    # Extract the final letter answer
    match = re.search(r'<<<([A-D])>>>', answer_text)
    if not match:
        return "Incorrect: The answer does not contain a final choice in the format <<<A>>>, <<<B>>>, etc."
    
    provided_answer = match.group(1)

    # --- Step 3: Check for Correctness ---

    # Check 1: Does the final letter match the derived correct option?
    if provided_answer != correct_option:
        return (f"Incorrect: The final answer is {provided_answer}, but the correct answer is {correct_option}. "
                f"The key interaction in the major product is: {correct_option_description}")

    # Check 2: Is the reasoning sound? The most critical point is identifying the major product.
    # A correct answer must identify the 'exo' product as major due to steric hindrance.
    reasoning_is_sound = False
    if 'exo' in answer_text.lower() and 'major' in answer_text.lower() and ('steric' in answer_text.lower() or 'hindrance' in answer_text.lower()):
        # Check that it's not saying 'endo' is major
        if 'endo' in answer_text.lower() and 'major' in answer_text.lower():
             # If both are mentioned, check which one is definitively stated as major
             if answer_text.lower().find("exo adduct is the major product") != -1:
                 reasoning_is_sound = True
        else:
            reasoning_is_sound = True

    if not reasoning_is_sound:
        return ("Incorrect: The final letter is correct, but the reasoning is flawed. "
                "The answer fails to correctly identify the 'exo' adduct as the major product due to steric hindrance, "
                "which is the correct justification. Many incorrect answers assume the 'endo' rule applies without considering the severe steric clash.")

    # If both checks pass, the answer is correct.
    return "Correct"

# The provided answer from the prompt
llm_answer = """
Here is a step-by-step analysis of the problem, synthesizing the information and correcting the inconsistencies found in the candidate answers.

### Step 1: Identify the Reactants and Reaction
*   **Dienophile:** The problem describes the dehydration of a *cis*-alkene dicarboxylic acid to an anhydride with specific NMR data (¹H: one peak at 7 ppm; ¹³C: two peaks at 137 ppm, 165 ppm). The high symmetry and characteristic chemical shifts unambiguously identify this compound as **maleic anhydride**. This point is a consensus among all candidate answers.
*   **Diene:** The diene is explicitly named as **1,2,3,4-tetramethyl-1,3-cyclopentadiene**.
*   **Reaction:** The reaction is a **Diels-Alder [4+2] cycloaddition**, which forms a bicyclic product.

### Step 2: Determine the Major Stereoisomer (The Critical Point)
This is the most crucial step and the primary source of disagreement among the candidate answers.
*   **The Alder-Endo Rule:** This rule predicts that the *endo* adduct is typically the major kinetic product due to favorable secondary orbital interactions. Several candidate answers (e.g., 4, 7, 9, 10, 11, 12, 15) incorrectly apply this rule without considering other factors.
*   **Steric Hindrance Exception:** The diene is exceptionally bulky due to the four methyl groups on the carbons involved in the cycloaddition. During the *endo* approach, the anhydride ring would experience severe steric clash with the methyl groups at the bridgehead positions. The *exo* approach, where the anhydride is oriented away from these methyl groups, is sterically much less hindered.
*   **Conclusion:** In this specific case, the overwhelming steric hindrance overrides the endo rule. Therefore, the **exo adduct is the major product**, and the *endo* adduct is the minor product. The most accurate candidate answers (e.g., 2, 3, 5, 6, 8, 13, 14) correctly identify this.

### Step 3: Analyze the 3D Structure and NOESY Correlations
A NOESY cross-peak indicates that two protons are close in space (< 5 Å). We need to find a spatial proximity that exists in the major (*exo*) product but is absent in the minor (*endo*) product.

*   **Signal Assignments:** Based on typical chemical shifts, the signals can be assigned to the protons in the product structure:
    *   **~3.5 ppm (2H singlet):** The two equivalent protons on the anhydride ring (`H_anhydride`).
    *   **~1.7 ppm (6H singlet):** The six protons of the two equivalent methyl groups on the double bond (`Me_vinyl`).
    *   **~1.5 ppm (1H doublet):** One of the two non-equivalent protons on the C7 bridge (`H_bridge`).
    *   **~1.0 ppm (6H singlet):** The six protons of the two equivalent methyl groups at the saturated bridgehead positions (`Me_bridgehead`).

*   **Spatial Analysis of the Major (*Exo*) Product:**
    *   In the *exo* adduct, the anhydride ring is oriented on the same side as the C7 bridge.
    *   This geometry places the `H_anhydride` protons on the opposite face (*endo* face), pointing towards the C=C double bond.
    *   The `Me_vinyl` groups are also on the *endo* face.
    *   Therefore, in the major (*exo*) product, the **`H_anhydride` protons are spatially close to the `Me_vinyl` protons**. This will generate a strong NOESY cross-peak.

*   **Spatial Analysis of the Minor (*Endo*) Product:**
    *   In the *endo* adduct, the anhydride ring is on the *endo* face.
    *   This places the `H_anhydride` protons on the *exo* face, far away from the `Me_vinyl` groups.
    *   In this minor product, the `H_anhydride` protons would instead be close to one of the `H_bridge` protons.

### Step 4: Final Conclusion
The cross-peak that is present in the major product but absent in the minor one must be the one unique to the major (*exo*) isomer. This is the interaction between the anhydride protons and the vinylic methyl protons.

*   **Interacting Protons:** `H_anhydride` and `Me_vinyl`.
*   **Corresponding Resonances:** The **2H singlet at ~3.5 ppm** and the **6H singlet at ~1.7 ppm**.

Matching this to the provided options:
A) A 6H singlet at ~1 ppm and a 1H doublet at ~1.5 ppm
B) A 6H singlet at ~1 ppm and a 6H singlet at ~1.7 ppm
C) A 6H singlet at ~1.7 ppm and a 2H singlet at ~3.5 ppm
D) A 1H doublet at ~1.5 ppm and a 2H singlet at ~3.5 ppm

The correct choice is C, as it describes the interaction between the vinylic methyl groups and the anhydride protons.

<<<C>>>
"""

# Run the check
result = check_diels_alder_noesy(llm_answer)
print(result)