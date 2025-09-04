import re

def check_diels_alder_noesy_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the answer to the Diels-Alder/NOESY problem.

    The function follows these logical steps:
    1. Defines the chemical entities and their properties (protons, signals).
    2. Determines the major product of the reaction by considering steric effects.
    3. Identifies the unique NOESY cross-peak expected for the major product.
    4. Maps this cross-peak to the given multiple-choice options.
    5. Compares the derived correct option with the provided answer.
    """

    # Step 1: Define chemical entities and map signals to options
    # These are the key proton groups in the product.
    # The 1H doublet at ~1.5 ppm is one of the C7 bridge protons. In the endo adduct,
    # the anhydride protons are close to the 'anti' bridge proton.
    proton_types = {
        "H_anhydride": "2H singlet at ~3.5 ppm",
        "Me_vinyl": "6H singlet at ~1.7 ppm",
        "Me_bridgehead": "6H singlet at ~1.0 ppm",
        "H_bridge": "1H doublet at ~1.5 ppm"
    }

    # Map the multiple-choice options to the interacting proton pairs
    options_to_protons = {
        'A': frozenset({"Me_bridgehead", "Me_vinyl"}),
        'B': frozenset({"Me_bridgehead", "H_bridge"}),
        'C': frozenset({"H_bridge", "H_anhydride"}),
        'D': frozenset({"Me_vinyl", "H_anhydride"})
    }

    # Step 2: Determine the major product and its characteristic NOESY interaction
    # The Alder-endo rule predicts the 'endo' product is major.
    # However, the diene (1,2,3,4-tetramethyl-1,3-cyclopentadiene) is extremely
    # sterically hindered. This steric clash in the 'endo' transition state
    # makes the 'exo' pathway more favorable.
    # Conclusion: The 'exo' adduct is the major product.
    major_product_type = "exo"

    # Define the key distinguishing spatial proximities for each isomer.
    # This is the interaction that is strong in one isomer and weak/absent in the other.
    distinguishing_proximities = {
        "endo": frozenset({"H_anhydride", "H_bridge"}),  # Corresponds to option C
        "exo": frozenset({"H_anhydride", "Me_vinyl"})    # Corresponds to option D
    }

    # Step 3: Identify the correct option based on the major product
    # The question asks for the cross-peak present in the major product.
    # Since the major product is 'exo', we look for the 'exo' proximity.
    correct_proton_pair = distinguishing_proximities[major_product_type]

    correct_option = None
    for option, protons in options_to_protons.items():
        if protons == correct_proton_pair:
            correct_option = option
            break

    if correct_option is None:
        # This case should not be reached with the current logic.
        return "Error in checker: Could not map the derived correct proton pair to an option."

    # Step 4: Extract the answer from the LLM's response and check it
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Could not parse the answer format from the text: '{llm_answer_text}'"
    
    llm_option = match.group(1)

    if llm_option == correct_option:
        return "Correct"
    else:
        reason = (
            f"The provided answer is '{llm_option}', but the correct answer is '{correct_option}'.\n\n"
            f"Reasoning:\n"
            f"1. The reaction is a Diels-Alder cycloaddition. The key decision is identifying the major product.\n"
            f"2. While the Alder-endo rule often applies, the diene in this problem (1,2,3,4-tetramethyl-1,3-cyclopentadiene) is exceptionally bulky. This severe steric hindrance favors the 'exo' transition state, making the 'exo' adduct the major product.\n"
            f"3. The question asks for a NOESY cross-peak present in the major product but absent in the minor. We must therefore analyze the 'exo' structure.\n"
            f"4. In the 'exo' adduct, the anhydride protons (~3.5 ppm) are spatially close to the vinylic methyl groups (~1.7 ppm). This proximity is absent in the 'endo' adduct.\n"
            f"5. This interaction corresponds to option {correct_option}. The provided answer '{llm_option}' is incorrect. An answer of 'C', for example, would incorrectly assume the 'endo' adduct is the major product."
        )
        return reason

# The final answer provided by the LLM to be checked
final_answer_from_llm = """
This problem requires a careful step-by-step analysis of a Diels-Alder reaction, including its stereochemical outcome and the interpretation of 2D NOESY NMR data.

### Step 1: Identify the Reactants

1.  **Dienophile (Anhydride):** The problem describes a *cis*-dicarboxylic acid that is dehydrated to an anhydride. The anhydride's NMR data is key:
    *   **¹H NMR:** A single peak at 7 ppm. This indicates all protons are chemically equivalent.
    *   **¹³C NMR:** Two peaks at 137 ppm (alkene carbons) and 165 ppm (carbonyl carbons). This indicates two sets of equivalent carbons.
    *   This high degree of symmetry, combined with the chemical shifts, unambiguously identifies the dienophile as **maleic anhydride**.

2.  **Diene:** The diene is explicitly named as **1,2,3,4-tetramethyl-1,3-cyclopentadiene**.

### Step 2: Analyze the Diels-Alder Reaction and its Products

The reaction is a [4+2] cycloaddition between the diene and dienophile, which produces a bicyclo[2.2.1]heptene derivative. Two stereoisomers are possible: the *endo* and *exo* adducts. The question states that a major and a minor product are formed.

### Step 3: Determine the Major Product (*Endo* vs. *Exo*)

*   **The Alder Endo Rule:** This rule states that the *endo* product is typically the major kinetic product due to stabilizing secondary orbital interactions.
*   **Steric Hindrance Exception:** This rule can be overridden by significant steric hindrance. The diene, 1,2,3,4-tetramethyl-1,3-cyclopentadiene, is extremely bulky. The four methyl groups create a crowded environment. The *endo* transition state would force the anhydride ring to approach from the same face as the two bridgehead methyl groups, leading to severe steric clash. The *exo* transition state, where the anhydride approaches from the less hindered face, is therefore favored.
*   **Conclusion:** In this specific case, due to severe steric hindrance, the **exo adduct is the major product**, and the *endo* adduct is the minor product.

### Step 4: Analyze the 3D Structure and Spatial Proximity for NOESY

NOESY spectroscopy detects protons that are close in space (< 5 Å). We need to find a pair of protons that are close in the major (*exo*) product but far apart in the minor (*endo*) product.

*   **Signal Assignments:**
    *   **~3.5 ppm (2H singlet):** The two equivalent protons on the anhydride ring (`H_anhydride`).
    *   **~1.7 ppm (6H singlet):** The six protons of the two equivalent methyl groups on the double bond (`Me_vinyl`).
    *   **~1.0 ppm (6H singlet):** The six protons of the two equivalent methyl groups at the bridgehead positions (`Me_bridgehead`).
    *   **~1.5 ppm (1H doublet):** One of the two non-equivalent protons on the C7 methylene bridge (`H_bridge`).

*   **Spatial Analysis:**
    *   **In the Major (*Exo*) Product:** The anhydride ring points away from the C=C double bond. The `H_anhydride` protons are on the *endo* face of the bicyclic system. This orientation places them very close to the `Me_vinyl` groups, which are also on the *endo* face.
    *   **In the Minor (*Endo*) Product:** The anhydride ring is tucked under the C=C double bond. The `H_anhydride` protons are on the *exo* face of the bicyclic system. This places them far away from the `Me_vinyl` groups.

### Step 5: Final Conclusion

The NOESY cross-peak that is present in the major product but absent in the minor product must arise from the unique spatial proximity in the *exo* isomer. This is the interaction between the anhydride protons and the vinylic methyl protons.

*   This corresponds to a cross-peak between the **2H singlet at ~3.5 ppm** (`H_anhydride`) and the **6H singlet at ~1.7 ppm** (`Me_vinyl`).

This matches option D.

<<<D>>>
"""

# Run the check
result = check_diels_alder_noesy_answer(final_answer_from_llm)
print(result)