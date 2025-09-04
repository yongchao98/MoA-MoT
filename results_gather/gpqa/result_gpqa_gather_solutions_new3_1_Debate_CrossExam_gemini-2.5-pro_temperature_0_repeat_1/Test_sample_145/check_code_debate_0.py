import re

def check_diels_alder_stereochemistry(question, options, final_answer_text):
    """
    Checks the correctness of an answer to a Diels-Alder stereochemistry problem.

    This function codifies the chemical principles for the reaction between
    5-fluorocyclopenta-1,3-diene and maleic anhydride.

    Args:
        question (str): The question text.
        options (dict): A dictionary mapping option letters to their IUPAC names.
        final_answer_text (str): The full text of the answer to be checked.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # --- Define Correct Chemical Principles and Nomenclature ---

    # 1. Endo/Exo Selectivity Principle
    # For this reaction, the endo adduct is the major kinetic product.
    # The endo skeleton has the (4S, 7R) or (4R, 7S) relative stereochemistry.
    # The exo skeleton has the (4R, 7R) or (4S, 7S) relative stereochemistry.
    # Let's map the options to endo/exo based on their stereodescriptors.
    # A: (3aR,4S,7R,7aS,...) -> 4S, 7R -> Endo
    # B: (3aR,4R,7S,7aS,...) -> 4R, 7S -> Exo
    # C: (3aR,4R,7S,7aS,...) -> 4R, 7S -> Exo
    # D: (3aR,4S,7R,7aS,...) -> 4S, 7R -> Endo
    endo_options = ['A', 'D']
    exo_options = ['B', 'C']
    correct_endo_exo_conclusion = f"The kinetically favored 'endo' rule applies, making the major product an endo adduct. This eliminates options {exo_options}."

    # 2. Facial Selectivity Principle
    # For a C5-substituent that is a first-row heteroatom (like F), electronic effects
    # favor syn-facial attack (dienophile approaches from the same side as the F).
    correct_facial_selectivity_principle = "Electronic effects from the C5-Fluorine substituent favor 'syn-facial attack'."

    # 3. Relationship between Attack Geometry and Adduct Geometry
    # This is a crucial and often confused point.
    # - syn-facial attack -> places the substituent and dienophile on OPPOSITE sides of the final ring system -> ANTI-adduct.
    # - anti-facial attack -> places the substituent and dienophile on the SAME side of the final ring system -> SYN-adduct.
    correct_attack_adduct_relation = "A 'syn-facial attack' results in an 'anti-adduct', where the fluorine and the anhydride ring are on opposite sides of the bicyclic system."

    # 4. Nomenclature of the Adduct
    # '8s' descriptor corresponds to the syn-adduct (F and anhydride on the same side).
    # '8r' descriptor corresponds to the anti-adduct (F and anhydride on opposite sides).
    # Option A has '8s' -> syn-adduct
    # Option D has '8r' -> anti-adduct
    correct_nomenclature_conclusion = "The 'anti-adduct' is assigned the '8r' stereodescriptor. The 'syn-adduct' is assigned '8s'."

    # 5. Final Correct Conclusion
    # Endo selectivity + syn-facial attack -> endo, anti-adduct -> Option D
    correct_final_product_name = options['D']
    correct_final_option = 'D'

    # --- Analyze the Provided Answer ---
    try:
        # Extract the reasoning and the final chosen option from the provided text
        reasoning = final_answer_text.split('<<<')[0]
        match = re.search(r'<<<([A-D])>>>', final_answer_text)
        if not match:
            return "Error: Could not find a final answer in the format <<<A>>>."
        chosen_option = match.group(1)
    except Exception as e:
        return f"Error parsing the provided answer: {e}"

    # Check the provided answer's logic against the correct principles
    errors = []

    # Check 1: Does the answer correctly identify the endo rule?
    if not ("endo" in reasoning.lower() and ("kinetic" in reasoning.lower() or "orbital" in reasoning.lower())):
        errors.append("The reasoning fails to correctly invoke the 'endo rule' as the basis for selectivity.")
    
    # Check 2: Does the answer correctly identify syn-facial attack?
    if "syn-facial attack" not in reasoning.lower():
        errors.append("The reasoning fails to identify that 'syn-facial attack' is electronically favored.")

    # Check 3: Does the answer correctly relate attack geometry to adduct geometry?
    # This is the most likely point of failure.
    if "syn-facial attack" in reasoning.lower() and "syn adduct" in reasoning.lower():
        error_msg = (
            "The answer makes a critical error in stereochemical reasoning. It incorrectly states that 'syn-facial attack' "
            "leads to a 'syn-adduct'. The correct relationship is that syn-facial attack (where the dienophile approaches "
            "from the same side as the substituent) results in an 'anti-adduct' (where the substituent and the anhydride "
            "ring are on opposite sides of the final product)."
        )
        errors.append(error_msg)

    # Check 4: Does the answer choose the correct final product?
    if chosen_option != correct_final_option:
        error_msg = (
            f"The answer incorrectly selects option {chosen_option}. The combination of endo selectivity and syn-facial attack "
            f"leads to the 'endo, anti-adduct'. This corresponds to option {correct_final_option} ({correct_final_product_name}), "
            f"not option {chosen_option}."
        )
        errors.append(error_msg)

    if errors:
        return "Incorrect. " + " ".join(errors)
    else:
        return "Correct"

# --- Execution ---
question_text = "5-fluorocyclopenta-1,3-diene is reacted with maleic anhydride. What is the major product?"
options_dict = {
    'A': "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
    'B': "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
    'C': "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
    'D': "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione"
}
# The provided answer from the prompt
llm_answer_text = """
Here is a step-by-step analysis to determine the major product of the reaction between 5-fluorocyclopenta-1,3-diene and maleic anhydride.

1.  **Identify the Reaction Type**: The reaction is a [4+2] cycloaddition, known as the Diels-Alder reaction. The reactants are a substituted cyclic diene (5-fluorocyclopenta-1,3-diene) and a reactive dienophile (maleic anhydride). This reaction forms a bridged bicyclic product.

2.  **Determine Endo/Exo Selectivity**: The Diels-Alder reaction is typically under kinetic control, and for cyclic dienes like cyclopentadiene reacting with dienophiles like maleic anhydride, the major product is the *endo* adduct. This preference arises from stabilizing secondary orbital interactions between the π-system of the diene and the electron-withdrawing carbonyl groups of the dienophile in the transition state. The *endo* product is formed faster, even though the *exo* product is often more thermodynamically stable.

3.  **Determine Facial Selectivity (Syn/Anti Attack)**: The fluorine atom at the C5 position of the diene makes the two faces of the diene non-equivalent. The dienophile can approach from the same side as the fluorine (*syn*-attack) or the opposite side (*anti*-attack).
    *   While steric hindrance would suggest that *anti*-attack is favored, it is a known principle in Diels-Alder reactions that electronegative substituents at the C5 position (like F, Cl, OAc) electronically favor a *syn*-facial attack. This "contrasteric" selectivity is due to stabilizing orbital interactions in the *syn* transition state.
    *   Therefore, the major pathway involves the dienophile approaching from the same side as the fluorine atom (*syn*-attack).

4.  **Predict the Major Product's Structure**: The major product is the one that forms fastest, which results from the combination of the two most favorable kinetic pathways: an **endo** transition state and **syn-facial attack**.
    *   The *endo* orientation places the anhydride ring "under" the diene system.
    *   The *syn*-attack means the dienophile approaches from the same side as the fluorine.
    *   The resulting product is the ***endo, syn*** **adduct**, where the fluorine atom (on the C8 bridge) and the anhydride ring are on the same side of the bicyclic system.

5.  **Match the Structure to the IUPAC Name**: The final step is to correlate the *endo, syn* structure with the correct IUPAC name among the choices.
    *   **Endo/Exo Ring Stereochemistry**: The stereochemical descriptors `(3aR,4S,7R,7aS)` represent the *endo* adduct skeleton. The descriptors `(3aR,4R,7S,7aS)` represent the *exo* skeleton. Since the *endo* product is major, we can eliminate options B and C.
    *   **Syn/Anti Fluorine Stereochemistry**: We are left with options A and D, which are both *endo* adducts but differ at the C8 position (`8s` vs. `8r`). The *syn* position places the fluorine atom on the same side as the anhydride ring. By assigning Cahn-Ingold-Prelog (CIP) priorities to the atoms around the C8 carbon in this *syn* arrangement, the configuration is determined to be **`s`**. The *anti* position would correspond to the `r` configuration.
    *   Therefore, the correct name for the major *endo, syn* product is **(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione**.

This corresponds to option A.

<<<A>>>
"""

# Run the check
result = check_diels_alder_stereochemistry(question_text, options_dict, llm_answer_text)
print(result)