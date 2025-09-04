import re

def check_diels_alder_noesy_answer(final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the Diels-Alder NOESY question.

    This function verifies the chain of logic in the provided answer against
    established chemical principles for this specific reaction.
    """
    # Step 1: Define the ground truth based on chemical principles
    # Options provided in the question
    options = {
        'A': {'6H singlet at ~1 ppm', '6H singlet at ~1.7 ppm'},
        'B': {'6H singlet at ~1.7 ppm', '2H singlet at ~3.5 ppm'},
        'C': {'1H doublet at ~1.5 ppm', '2H singlet at ~3.5 ppm'},
        'D': {'6H singlet at ~1 ppm', '1H doublet at ~1.5 ppm'}
    }

    # Proton signal assignments to chemical groups
    proton_assignments = {
        'anhydride_protons': '2H singlet at ~3.5 ppm',
        'vinylic_methyls': '6H singlet at ~1.7 ppm',
        'bridgehead_methyls': '6H singlet at ~1.0 ppm',
        'bridge_proton': '1H doublet at ~1.5 ppm'
    }

    # Spatial proximities to anhydride protons in each isomer, determining the unique NOE.
    # In the exo adduct, anhydride protons are on the endo face, near the vinylic methyls.
    # In the endo adduct, anhydride protons are on the exo face, near the bridge proton.
    noe_proximities = {
        'exo': ['vinylic_methyls'],
        'endo': ['bridge_proton']
    }

    # Step 2: Parse the provided answer to extract the choice and reasoning
    try:
        chosen_option_match = re.search(r'<<<([A-D])>>>', final_answer_text)
        if not chosen_option_match:
            return "Failure: Could not find a final answer in the format <<<A>>>."
        chosen_option = chosen_option_match.group(1)
        reasoning_text = final_answer_text.lower()
    except Exception as e:
        return f"Failure: Error parsing the provided answer: {e}"

    # Step 3: Evaluate the reasoning in the provided answer
    
    # Check 3a: Does the reasoning identify the major product?
    # The most chemically sound argument for this specific reaction is that the exo product is major due to sterics.
    is_exo_major = 'exo' in reasoning_text and 'major' in reasoning_text
    is_endo_major = 'endo' in reasoning_text and 'major' in reasoning_text

    if is_exo_major:
        major_isomer = 'exo'
    elif is_endo_major:
        major_isomer = 'endo'
    else:
        return "Incorrect: The reasoning fails to clearly state whether the 'endo' or 'exo' product is the major one, which is a critical first step."

    # Check 3b: Based on the identified major product, does the reasoning correctly identify the NOE interaction?
    expected_interacting_groups = noe_proximities[major_isomer]
    
    mentions_vinylic_interaction = 'vinylic methyl' in reasoning_text and ('anhydride proton' in reasoning_text or 'anhydride ring' in reasoning_text)
    mentions_bridge_interaction = 'bridge proton' in reasoning_text and ('anhydride proton' in reasoning_text or 'anhydride ring' in reasoning_text)

    if major_isomer == 'exo':
        if not mentions_vinylic_interaction:
            return "Incorrect: The reasoning correctly identifies 'exo' as major but fails to deduce the correct spatial proximity between anhydride protons and vinylic methyls."
        expected_signals = {proton_assignments['anhydride_protons'], proton_assignments['vinylic_methyls']}
        correct_option_for_reasoning = 'B'

    elif major_isomer == 'endo':
        # The correct interaction for endo is with the bridge proton.
        if mentions_vinylic_interaction:
            return "Incorrect: The reasoning claims 'endo' is major but incorrectly states an interaction with vinylic methyls. In the endo isomer, these groups are far apart."
        if not mentions_bridge_interaction:
            return "Incorrect: The reasoning identifies 'endo' as major but fails to deduce the correct spatial proximity between anhydride protons and the bridge proton."
        expected_signals = {proton_assignments['anhydride_protons'], proton_assignments['bridge_proton']}
        correct_option_for_reasoning = 'C'

    # Step 4: Check if the chosen option matches the logical conclusion from the reasoning
    def normalize_signal(s):
        return s.replace(' at ', '').replace('~', '').strip()

    normalized_expected_signals = {normalize_signal(s) for s in expected_signals}
    normalized_chosen_signals = {normalize_signal(s) for s in options[chosen_option]}

    if normalized_expected_signals == normalized_chosen_signals:
        if chosen_option == correct_option_for_reasoning:
            return "Correct"
        else:
            # This case should not be reached if logic is sound, but as a safeguard:
            return f"Incorrect: Internal inconsistency. The reasoning leads to signals matching option {correct_option_for_reasoning}, but option {chosen_option} was selected."
    else:
        return (f"Incorrect: The final choice is inconsistent with the reasoning. The reasoning for a major '{major_isomer}' isomer "
                f"points to an interaction corresponding to option {correct_option_for_reasoning}. "
                f"However, option {chosen_option} was selected.")

# The final answer provided by the user to check
final_answer_to_check = """
### Step-by-step Analysis

1.  **Identify the Reactants:**
    *   **Dienophile:** The problem describes a *cis*-dicarboxylic acid that is dehydrated to an anhydride. The NMR data for this anhydride (¹H: one peak at 7 ppm; ¹³C: two peaks at 137 ppm and 165 ppm) indicates a highly symmetrical molecule. This data is a perfect match for **maleic anhydride**, which is formed from the dehydration of maleic acid. The two vinylic protons are equivalent, giving one ¹H signal. The two vinylic carbons are equivalent, and the two carbonyl carbons are equivalent, giving two ¹³C signals.
    *   **Diene:** The diene is explicitly named as **1,2,3,4-tetramethyl-1,3-cyclopentadiene**.

2.  **Identify the Reaction and Products:**
    *   The reaction is a **Diels-Alder reaction**, a [4+2] cycloaddition between the diene and the dienophile.
    *   This reaction forms a bicyclo[2.2.1]heptene derivative.
    *   The formation of a major and a minor product indicates the formation of two stereoisomers: the ***endo*** and ***exo*** adducts.

3.  **Determine the Major Product and its Stereochemistry:**
    *   The **Alder endo rule** states that the *endo* product is typically the major, kinetically favored product due to stabilizing secondary orbital interactions.
    *   However, this rule can be overridden by severe steric hindrance. The diene, 1,2,3,4-tetramethyl-1,3-cyclopentadiene, is extremely bulky. The four methyl groups create significant steric clash in the *endo* transition state. The *exo* transition state, where the anhydride ring approaches from the side away from the C7 methylene bridge, is sterically less hindered.
    *   Therefore, in this specific case, the **major product is the *exo* adduct**, and the minor product is the *endo* adduct.
    *   Let's define the stereochemistry based on the standard convention for bicyclo[2.2.1] systems:
        *   ***Exo* Adduct (Major Product):** The anhydride ring is positioned *anti* (on the opposite side) to the C7 bridge. In this orientation, the protons on the anhydride portion (H₅, H₆) are on the *endo* face of the bicyclic system, tucked "under" the C=C double bond.
        *   ***Endo* Adduct (Minor Product):** The anhydride ring is positioned *syn* (on the same side) to the C7 bridge. The anhydride protons are on the *exo* face, pointing away from the C=C double bond.

4.  **Analyze the NOESY Spectrum:**
    *   NOESY (Nuclear Overhauser Effect Spectroscopy) detects protons that are close in space (< 5 Å). A cross-peak indicates spatial proximity.
    *   The question states a cross-peak is **present in the major product** but absent in the minor product.
    *   Based on our analysis in Step 3:
        *   In the **major (*exo*) product**, the anhydride protons are on the *endo* face, which brings them very close to the vinylic methyl groups (at C2 and C3), which are also on the *endo* face. A strong NOESY cross-peak is expected between them.
        *   In the **minor (*endo*) product**, the anhydride protons are on the *exo* face, far away from the vinylic methyl groups. This cross-peak would be absent or very weak.
    *   Therefore, the distinguishing cross-peak connects the anhydride protons and the vinylic methyl protons.

5.  **Assign NMR Signals and Select the Correct Option:**
    *   Based on typical chemical shifts for this type of structure:
        *   **~1.0 ppm (6H singlet):** The two equivalent methyl groups at the saturated bridgehead positions.
        *   **~1.7 ppm (6H singlet):** The two equivalent methyl groups on the C=C double bond (vinylic methyls).
        *   **~3.5 ppm (2H singlet):** The two equivalent protons on the anhydride ring, deshielded by the adjacent carbonyl groups.
        *   **~1.5 ppm (1H doublet):** One of the two non-equivalent protons on the C7 bridge.
    *   The cross-peak in the major product connects the vinylic methyl protons and the anhydride protons. This corresponds to the signals at **~1.7 ppm (6H singlet)** and **~3.5 ppm (2H singlet)**.

6.  **Final Check:**
    *   The option that describes a cross-peak between a 6H singlet at ~1.7 ppm and a 2H singlet at ~3.5 ppm is option B.

<<<B>>>
"""

# Run the check
result = check_diels_alder_noesy_answer(final_answer_to_check)
print(result)