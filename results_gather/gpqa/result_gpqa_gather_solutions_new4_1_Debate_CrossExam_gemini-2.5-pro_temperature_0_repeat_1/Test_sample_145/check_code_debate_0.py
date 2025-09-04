import re

def check_answer(selected_answer_letter: str, detailed_reasoning: str) -> str:
    """
    Checks the correctness of the answer for the Diels-Alder reaction.

    The function verifies the answer based on two key stereochemical principles:
    1. The Alder Endo Rule: The kinetic product is the 'endo' adduct.
    2. Facial Selectivity for 5-F-cyclopentadiene: Electronic effects favor 'syn'-facial
       attack, which results in an 'anti' product (fluorine is anti to the anhydride ring).

    Args:
        selected_answer_letter: The letter of the chosen answer (e.g., 'B').
        detailed_reasoning: The reasoning provided for the final answer.

    Returns:
        "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # Step 1: Define the properties of each answer option based on IUPAC nomenclature.
    # 'endo_exo': Describes the orientation of the anhydride ring.
    # 'syn_anti': Describes the position of the fluorine relative to the anhydride ring.
    # '8s' -> syn (same side), '8r' -> anti (remote/opposite side)
    # '(4S,7R)' -> endo, '(4R,7S)' -> exo
    options_properties = {
        'A': {'endo_exo': 'endo', 'syn_anti': 'syn'},
        'B': {'endo_exo': 'endo', 'syn_anti': 'anti'},
        'C': {'endo_exo': 'exo',  'syn_anti': 'syn'},
        'D': {'endo_exo': 'exo',  'syn_anti': 'anti'},
    }

    # Step 2: Define the expected properties of the major product based on chemical principles.
    # Principle 1: The Alder Endo Rule predicts the 'endo' adduct is the major kinetic product.
    expected_endo_exo = 'endo'
    # Principle 2: For a C5-Fluorine substituent, electronic effects favor 'syn'-facial attack,
    # which leads to an 'anti' configuration in the final product.
    expected_syn_anti = 'anti'

    # Step 3: Check if the selected answer exists in our defined options.
    if selected_answer_letter not in options_properties:
        return f"Invalid answer option '{selected_answer_letter}'. Please choose from A, B, C, or D."

    # Step 4: Retrieve the properties of the selected answer.
    selected_properties = options_properties[selected_answer_letter]

    # Step 5: Verify the selected answer against the chemical principles.
    # Check Endo/Exo selectivity
    if selected_properties['endo_exo'] != expected_endo_exo:
        return (f"Incorrect. The answer '{selected_answer_letter}' corresponds to an '{selected_properties['endo_exo']}' adduct. "
                f"The major product of this Diels-Alder reaction should be the '{expected_endo_exo}' adduct "
                f"due to the Alder Endo Rule (kinetic control).")

    # Check Syn/Anti selectivity
    if selected_properties['syn_anti'] != expected_syn_anti:
        return (f"Incorrect. The answer '{selected_answer_letter}' corresponds to a '{selected_properties['syn_anti']}' product. "
                f"The major product should be the '{expected_syn_anti}' adduct. This results from the favored "
                f"'syn'-facial attack of the dienophile on the 5-fluorocyclopentadiene, placing the fluorine "
                f"and anhydride ring on opposite sides.")

    # Step 6: If all checks pass, the answer is correct.
    return "Correct"

# --- Execution ---
# Extract the final answer and reasoning from the provided text.
final_answer_text = """
Here is a step-by-step analysis of the problem, considering the points of contention and potential errors highlighted in the candidate answers.

### Step 1: Identify the Reaction and Key Principles

The reaction is a Diels-Alder cycloaddition between 5-fluorocyclopenta-1,3-diene and maleic anhydride. The major product is the one formed fastest, i.e., the kinetic product. Its stereochemistry is determined by two main factors:
1.  **Endo/Exo Selectivity:** The orientation of the dienophile (maleic anhydride) relative to the diene.
2.  **Syn/Anti Facial Selectivity:** The face of the diene that the dienophile attacks, relative to the fluorine substituent.

### Step 2: Analyze Endo/Exo Selectivity

The Diels-Alder reaction between cyclopentadiene and maleic anhydride is a classic textbook example of the **Alder Endo Rule**. This rule states that the kinetic product is the *endo* adduct, where the electron-withdrawing groups of the dienophile are tucked "under" the diene in the transition state. This is due to stabilizing secondary orbital interactions. While a C5 substituent can sometimes alter the selectivity, the *endo* pathway is generally still favored under kinetic control.

*   **Identifying the *Endo* Adducts:** The options are presented in two pairs of diastereomers that differ in their core stereochemistry.
    *   Options A and B share the core: `(3aR,4S,7R,7aS)`
    *   Options C and D share the core: `(3aR,4R,7S,7aS)`
    One of these patterns describes the *endo* adduct, and the other describes the *exo* adduct. A careful assignment of Cahn-Ingold-Prelog priorities confirms that the `(3aR,4S,7R,7aS)` skeleton corresponds to the **endo** adduct. Therefore, the major product must be either **A** or **B**.

### Step 3: Analyze Syn/Anti Facial Selectivity

This is the most critical and often misunderstood aspect. The dienophile can approach the diene from the same face as the C5-fluorine (*syn*-attack) or the opposite face (*anti*-attack).

*   **Steric Argument (Incorrect for Fluorine):** A simple steric model would suggest the dienophile should avoid the fluorine atom and perform an *anti*-attack.
*   **Electronic Argument (Correct for Fluorine):** For small, electronegative substituents at the C5 position (like F or OMe), electronic effects outweigh sterics. It is well-established in physical organic chemistry that such substituents favor **syn-facial attack**. This "contrasteric" selectivity is due to stabilizing orbital interactions in the transition state.

Therefore, the major reaction pathway involves the dienophile approaching the diene from the **same side as the fluorine atom**.

### Step 4: Determine the Final Product Structure

We must combine the outcomes of the two selective steps: an *endo* transition state and a *syn*-facial attack.

1.  Visualize the diene with the C5-fluorine pointing "up".
2.  The dienophile approaches from "up" as well (*syn*-attack).
3.  For the transition state to be *endo*, the anhydride's carbonyl groups must point towards the diene's double bonds. This orients the bulk of the anhydride ring "down", away from the diene.
4.  As the new bonds form, the final product has the C8 bridge (with the fluorine) pointing "up" and the anhydride ring fused in the "down" position.

This means the fluorine atom and the anhydride ring are on **opposite sides** of the final bicyclic system. This product is called the ***anti-endo*** **adduct**.

### Step 5: Match the Structure to the IUPAC Name

We have determined the major product is the *anti-endo* adduct and that it must be either option A or B. The final step is to determine the descriptor for the fluorine at position 8.

*   **A) (...,8s)-...**
*   **B) (...,8r)-...**

In the nomenclature for these bridged systems, the descriptors `s` and `r` are used to denote the relative position of a substituent on the C8 bridge.
*   `8s` (from *syn*) indicates the substituent is on the same side as the principal ring system (the anhydride ring).
*   `8r` (from *remote* or related to *anti*) indicates the substituent is on the side opposite the principal ring system.

Since our major product is the *anti-endo* adduct, the fluorine is *anti* (opposite) to the anhydride ring. Therefore, the correct descriptor is **8r**.

Combining all points:
*   Selectivity: **Endo** -> `(3aR,4S,7R,7aS)` skeleton.
*   Selectivity: **Anti** (relative position in product) -> `8r` descriptor.

The correct name for the major product is **(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione**.

<<<B>>>
"""

match = re.search(r'<<<([A-D])>>>', final_answer_text)
if match:
    answer_letter = match.group(1)
    result = check_answer(answer_letter, final_answer_text)
    print(result)
else:
    print("Could not find a final answer in the format <<<A>>>.")
