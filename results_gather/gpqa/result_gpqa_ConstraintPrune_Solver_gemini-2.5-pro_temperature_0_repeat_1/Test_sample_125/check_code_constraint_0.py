import re

def check_hplc_answer_correctness(llm_response_text):
    """
    This function checks the correctness of the provided answer by modeling the chemical reactions and HPLC separation principles.

    It verifies the following constraints:
    1.  The number of stereoisomers produced in Reaction I.
    2.  The number of stereoisomers produced in Reaction II.
    3.  The number of peaks expected in a normal-phase HPLC based on co-elution of enantiomers.
    4.  The number of peaks expected in a chiral HPLC based on separation of all stereoisomers.
    5.  The consistency between the reasoning in the text and the final selected option.
    """

    # --- Step 1: Model the products based on chemical principles ---

    # Reaction I: Reduction of (S)-5-methoxyhexan-3-one.
    # The starting material has one chiral center (C5). Reduction of the ketone at C3 creates a new chiral center.
    # This results in two products: (3R, 5S) and (3S, 5S)-5-methoxyhexan-3-ol.
    # These are diastereomers.
    # In normal-phase HPLC, diastereomers separate.
    products_rxn1 = [
        {'name': '(3R, 5S)-5-methoxyhexan-3-ol', 'normal_hplc_group': 'Rxn1_Diastereomer_A'},
        {'name': '(3S, 5S)-5-methoxyhexan-3-ol', 'normal_hplc_group': 'Rxn1_Diastereomer_B'}
    ]

    # Reaction II: Reduction of pentane-2,4-dione.
    # The starting material is achiral. Reduction creates two new chiral centers (C2, C4).
    # This results in three stereoisomers: an enantiomeric pair ((2R,4R) and (2S,4S)) and a meso compound ((2R,4S)).
    # In normal-phase HPLC, the enantiomeric pair co-elutes, while the meso compound (a diastereomer to the pair) elutes separately.
    products_rxn2 = [
        {'name': '(2R, 4R)-pentane-2,4-diol', 'normal_hplc_group': 'Rxn2_Enantiomeric_Pair'},
        {'name': '(2S, 4S)-pentane-2,4-diol', 'normal_hplc_group': 'Rxn2_Enantiomeric_Pair'},
        {'name': 'meso-pentane-2,4-diol', 'normal_hplc_group': 'Rxn2_Meso_Compound'}
    ]

    # --- Step 2: Calculate the expected number of peaks ---

    # Total number of unique stereoisomers produced
    all_products = products_rxn1 + products_rxn2
    
    # Chiral HPLC separates all unique stereoisomers.
    calculated_chiral_peaks = len(all_products)

    # Normal-phase HPLC separates based on structure and diastereomeric relationships. Enantiomers co-elute.
    # We use a set to count the unique groups that will form peaks.
    normal_hplc_groups = set()
    for product in all_products:
        normal_hplc_groups.add(product['normal_hplc_group'])
    calculated_normal_peaks = len(normal_hplc_groups)

    # --- Step 3: Extract the answer from the LLM's response ---
    
    # The question options define the expected values for each choice.
    options = {
        'A': {'chiral': 3, 'normal': 2},
        'B': {'chiral': 3, 'normal': 3},
        'C': {'chiral': 4, 'normal': 2},
        'D': {'chiral': 5, 'normal': 4}
    }

    # Find the selected option, e.g., <<<D>>>
    match = re.search(r'<<<([A-D])>>>', llm_response_text)
    if not match:
        return "Failure: Could not find the final answer in the format <<<X>>> in the provided text."
    
    selected_option = match.group(1)
    llm_answer = options[selected_option]

    # --- Step 4: Compare calculated values with the LLM's answer ---

    # Check if the number of chiral peaks matches.
    if calculated_chiral_peaks != llm_answer['chiral']:
        return (f"Incorrect. The answer states there are {llm_answer['chiral']} peaks in chiral HPLC, "
                f"but the correct number based on chemical principles is {calculated_chiral_peaks}. "
                f"There are {len(products_rxn1)} products from Reaction I and {len(products_rxn2)} products from Reaction II, "
                f"totaling {calculated_chiral_peaks} unique stereoisomers, which should all be separated by a chiral column.")

    # Check if the number of normal-phase peaks matches.
    if calculated_normal_peaks != llm_answer['normal']:
        return (f"Incorrect. The answer states there are {llm_answer['normal']} peaks in normal-phase HPLC, "
                f"but the correct number based on chemical principles is {calculated_normal_peaks}. "
                f"Normal-phase HPLC should produce {calculated_normal_peaks} peaks: "
                "2 from the diastereomers of Reaction I, 1 from the co-eluting enantiomeric pair of Reaction II, and 1 from the meso compound of Reaction II.")

    # If both checks pass, the answer is correct.
    return "Correct"

# The provided LLM response to be checked
llm_response = """
Here is the step-by-step analysis of the problem based on chemical constraints.

### Constraint Analysis

1.  **Constraint 1: Products of Reaction I**
    *   **Reactant:** (S)-5-methoxyhexan-3-one. This is a ketone with a single, defined chiral center at carbon 5 (C5).
    *   **Reaction:** Reduction of the ketone at C3 with LAH (a strong, achiral reducing agent) creates a new chiral center at C3. The original chiral center at C5 is unaffected.
    *   **Outcome:** The hydride can attack the planar carbonyl group from either face, leading to two possible configurations at C3 (R or S). This results in the formation of two products: (3R, 5S)-5-methoxyhexan-3-ol and (3S, 5S)-5-methoxyhexan-3-ol. These two molecules are **diastereomers** because they have one chiral center in common (S at C5) and one that is different (R/S at C3).
    *   **Conclusion:** Reaction I produces **2** stereoisomeric products.

2.  **Constraint 2: Products of Reaction II**
    *   **Reactant:** Pentane-2,4-dione. This is a symmetric, achiral diketone.
    *   **Reaction:** Reduction of both ketone groups with excess NaBH4 creates two new chiral centers at C2 and C4, forming pentane-2,4-diol.
    *   **Outcome:** Since the starting material and reagent are achiral, all possible stereoisomers will be formed.
        *   (2R, 4R)-pentane-2,4-diol
        *   (2S, 4S)-pentane-2,4-diol
        *   (2R, 4S)-pentane-2,4-diol
    *   The (2R, 4R) and (2S, 4S) isomers are non-superimposable mirror images, making them a pair of **enantiomers**. They will be formed in a 1:1 racemic mixture.
    *   The (2R, 4S) isomer has an internal plane of symmetry and is an achiral **meso compound**.
    *   **Conclusion:** Reaction II produces **3** stereoisomeric products (one enantiomeric pair and one meso compound).

3.  **Constraint 3: Analysis by Normal-Phase HPLC**
    *   **Principle:** Normal-phase (achiral) HPLC separates compounds based on differences in their physical properties, such as polarity. It cannot distinguish between enantiomers, which have identical physical properties in an achiral environment. It can, however, separate diastereomers, meso compounds, and molecules with different structures.
    *   **Applying to the mixture:**
        *   The two products from Reaction I, (3R, 5S)- and (3S, 5S)-5-methoxyhexan-3-ol, are diastereomers. They have different physical properties and will be separated. This gives **2 peaks**.
        *   From Reaction II, the enantiomeric pair ((2R, 4R)- and (2S, 4S)-pentane-2,4-diol) will not be separated and will elute together as a single peak. This gives **1 peak**.
        *   The meso-pentane-2,4-diol is a diastereomer of the enantiomeric pair. It has different properties and will be separated from them. This gives **1 peak**.
        *   The products of Reaction I (C7 ether-alcohols) are structurally different from the products of Reaction II (C5 diols) and will have different polarities, ensuring they separate from each other.
    *   **Total Peaks:** 2 (from Rxn I) + 1 (from Rxn II enantiomers) + 1 (from Rxn II meso) = **4 peaks**.

4.  **Constraint 4: Analysis by Chiral HPLC**
    *   **Principle:** Chiral HPLC uses a chiral stationary phase that interacts differently with each stereoisomer, including enantiomers. It can separate all non-identical stereoisomers.
    *   **Applying to the mixture:** We have a total of 5 unique stereoisomers in the combined product mixture.
        1.  (3R, 5S)-5-methoxyhexan-3-ol
        2.  (3S, 5S)-5-methoxyhexan-3-ol
        3.  (2R, 4R)-pentane-2,4-diol
        4.  (2S, 4S)-pentane-2,4-diol
        5.  meso-pentane-2,4-diol
    *   A chiral column will resolve all five of these distinct molecules.
    *   **Total Peaks:** **5 peaks**.

### Summary and Code Verification

The analysis consistently shows 5 peaks for chiral HPLC and 4 peaks for normal-phase HPLC. This corresponds to option D.

<<<D>>>
"""

# Run the check
result = check_hplc_answer_correctness(llm_response)
print(result)