import re

def check_correctness_of_chemistry_answer(question_text, final_answer_text):
    """
    Checks the correctness of the answer to the organic chemistry HPLC question.

    This function programmatically models the chemical reactions and chromatographic principles
    to determine the correct number of peaks and compares it with the provided answer.
    """
    # --- Step 1: Model the chemical reactions and their stereochemical outcomes ---

    # Reaction I: Reduction of (S)-5-methoxyhexan-3-one
    # A chiral starting material with one stereocenter (C5) is reacted to create a new stereocenter (C3).
    # This creates a pair of diastereomers.
    # - Number of unique stereoisomers = 2
    # - Peaks in normal HPLC (diastereomers separate) = 2
    # - Peaks in chiral HPLC (diastereomers separate) = 2
    rxn1_normal_peaks = 2
    rxn1_chiral_peaks = 2

    # Reaction II: Reduction of pentane-2,4-dione
    # An achiral starting material is reacted to create two new stereocenters (C2, C4).
    # This results in a pair of enantiomers and a meso compound.
    # - Number of unique stereoisomers = 3 (R,R enantiomer; S,S enantiomer; R,S meso)
    # - Peaks in normal HPLC (enantiomers co-elute, meso separates) = 1 + 1 = 2
    # - Peaks in chiral HPLC (all 3 stereoisomers separate) = 3
    rxn2_normal_peaks = 2
    rxn2_chiral_peaks = 3

    # --- Step 2: Calculate the total expected peaks for the combined mixture ---

    # The products from Reaction I and Reaction II are structurally different (different molecular formulas)
    # and will not co-elute in either HPLC setup.
    expected_normal_peaks = rxn1_normal_peaks + rxn2_normal_peaks
    expected_chiral_peaks = rxn1_chiral_peaks + rxn2_chiral_peaks

    # --- Step 3: Parse the options from the question and the provided answer ---

    # Regex to find options like "A) X Peaks in chiral HPLC and Y peaks in normal-phase HPLC"
    # It handles variations in wording and the order of chiral/normal.
    option_pattern = re.compile(
        r'([A-D])\)\s*(?P<val1>\d+)\s*Peaks?\s*in\s*(?P<type1>chiral|normal-phase)\s*HPLC\s*and\s*(?P<val2>\d+)\s*peaks?\s*in\s*(?P<type2>normal-phase|chiral)\s*HPLC',
        re.IGNORECASE
    )
    # Regex for options like "X peaks in both..."
    both_pattern = re.compile(
        r'([A-D])\)\s*(?P<val>\d+)\s*peaks?\s*in\s*both\s*chiral\s*HPLC\s*and\s*normal-phase\s*HPLC',
        re.IGNORECASE
    )

    options = {}
    for line in question_text.split('\n'):
        line = line.strip()
        match = option_pattern.match(line)
        if match:
            key = match.group(1)
            data = {}
            data[match.group('type1').lower().replace('-phase', '')] = int(match.group('val1'))
            data[match.group('type2').lower().replace('-phase', '')] = int(match.group('val2'))
            options[key] = data
        else:
            match_both = both_pattern.match(line)
            if match_both:
                key = match_both.group(1)
                val = int(match_both.group('val'))
                options[key] = {'chiral': val, 'normal': val}

    # Extract the final letter answer from the provided text
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Failure: Could not parse the final answer in the format <<<X>>>."
    llm_answer_key = match.group(1)

    if llm_answer_key not in options:
        return f"Failure: The final answer '{llm_answer_key}' is not a valid option key (A, B, C, D)."

    llm_selected_option = options[llm_answer_key]
    llm_chiral_peaks = llm_selected_option.get('chiral')
    llm_normal_peaks = llm_selected_option.get('normal')

    # --- Step 4: Compare calculated results with the provided answer ---

    is_correct = (llm_chiral_peaks == expected_chiral_peaks) and \
                 (llm_normal_peaks == expected_normal_peaks)

    if is_correct:
        return "Correct"
    else:
        # Find the correct option key for the error message
        correct_key = "None"
        for key, value in options.items():
            if value.get('chiral') == expected_chiral_peaks and value.get('normal') == expected_normal_peaks:
                correct_key = key
                break

        reason = (
            f"The provided answer '{llm_answer_key}' is incorrect.\n"
            f"Constraint Analysis:\n"
            f"1. Reaction I produces 2 diastereomers.\n"
            f"2. Reaction II produces 3 unique stereoisomers (1 enantiomeric pair + 1 meso compound).\n"
            f"3. Normal-phase HPLC separates diastereomers but not enantiomers. Expected peaks = 2 (from Rxn I) + 2 (from Rxn II) = {expected_normal_peaks}.\n"
            f"4. Chiral HPLC separates all unique stereoisomers. Expected peaks = 2 (from Rxn I) + 3 (from Rxn II) = {expected_chiral_peaks}.\n\n"
            f"Conclusion:\n"
            f"The correct answer requires {expected_chiral_peaks} chiral peaks and {expected_normal_peaks} normal-phase peaks, which corresponds to option '{correct_key}'.\n"
            f"The selected answer '{llm_answer_key}' claims {llm_chiral_peaks} chiral peaks and {llm_normal_peaks} normal-phase peaks, which violates the constraints."
        )
        return reason

# --- Execution ---
# Extract the relevant parts from the user's prompt
question_full_text = """
An organic chemist performs two reactions.
Reaction I:
(S)-5-methoxyhexan-3-one is treated with LAH followed by an acidic workup.
Reaction II:
Pentane-2,4-dione is treated with excess NaBH4 followed by an acidic workup.
The progress of the reactions is monitored using TLC. After the reactions are 100% completed, the products from both reactions are combined, and run on both a normal-phase HPLC column and an HPLC column packed with the chiral stationary phase.
How many peaks does the chemist observe in the chromatograms?

A) 4 Peaks in chiral HPLC and 2 peaks in normal-phase HPLC
B) 5 Peaks in chiral HPLC and 4 peaks in normal-phase HPLC
C) 3 Peaks in chiral HPLC and 2 peaks in normal-phase HPLC
D) 3 peaks in both chiral HPLC and normal-phase HPLC
"""

final_answer_from_llm = """
Here is a step-by-step analysis to determine the final answer.

### **Step 1: Analyze the Products of Reaction I**

*   **Starting Material:** (S)-5-methoxyhexan-3-one. This is a chiral molecule with a single, pre-existing stereocenter at the C5 position, which has an (S) configuration.
*   **Reaction:** The ketone group at the C3 position is reduced to a secondary alcohol using Lithium Aluminum Hydride (LAH), an achiral reducing agent. This reduction creates a new stereocenter at the C3 position.
*   **Products:** The original stereocenter at C5 is not involved in the reaction and its (S) configuration is retained. The reduction of the planar carbonyl group at C3 can occur from two different faces, leading to both (R) and (S) configurations at the new C3 stereocenter. This results in two products:
    1.  (3R, 5S)-5-methoxyhexan-3-ol
    2.  (3S, 5S)-5-methoxyhexan-3-ol
*   **Relationship:** These two products are **diastereomers** because they have two stereocenters, one of which is the same (5S) and one of which is different (3R vs. 3S).
*   **Conclusion for Reaction I:** Two distinct diastereomeric products are formed.

### **Step 2: Analyze the Products of Reaction II**

*   **Starting Material:** Pentane-2,4-dione. This is a symmetric and **achiral** molecule.
*   **Reaction:** Both ketone groups are reduced to secondary alcohols using excess Sodium Borohydride (NaBH₄), an achiral reducing agent. This creates two new stereocenters at the C2 and C4 positions.
*   **Products:** Since the starting material and reagent are achiral, all possible stereoisomers are formed. The product is pentane-2,4-diol.
    1.  **(2R, 4R)-pentane-2,4-diol** and **(2S, 4S)-pentane-2,4-diol**: These two are non-superimposable mirror images, forming a pair of **enantiomers**. They are produced in a 1:1 racemic mixture.
    2.  **(2R, 4S)-pentane-2,4-diol**: This molecule has an internal plane of symmetry, making it an achiral **meso compound**.
*   **Relationship:** The meso compound is a diastereomer of the enantiomeric pair.
*   **Conclusion for Reaction II:** Three distinct stereoisomers are formed (a pair of enantiomers and a meso compound).

### **Step 3: Analyze the Combined Mixture with Normal-Phase HPLC**

*   **Principle:** A normal-phase HPLC column is achiral. It separates compounds based on differences in physical properties like polarity. It can separate constitutional isomers and diastereomers, but it **cannot** separate enantiomers.
*   **Peak Count:**
    *   **From Reaction I:** The two products are diastereomers. Diastereomers have different physical properties and will be separated. This gives **2 peaks**.
    *   **From Reaction II:** The products are an enantiomeric pair and a meso compound. The enantiomeric pair has identical physical properties in an achiral environment and will elute together as **1 peak**. The meso compound is a diastereomer of the pair, has different physical properties, and will elute as a separate, **second peak**. This gives a total of **2 peaks** for the products of Reaction II.
    *   **Overall:** The products from Reaction I (C7 ether-alcohols) are structurally different from those of Reaction II (C5 diols) and will not co-elute. Therefore, the total number of peaks is 2 (from Reaction I) + 2 (from Reaction II) = **4 peaks**.

### **Step 4: Analyze the Combined Mixture with Chiral HPLC**

*   **Principle:** A chiral HPLC column has a chiral stationary phase that can interact differently with all stereoisomers, allowing for the separation of **both diastereomers and enantiomers**.
*   **Peak Count:**
    *   The number of peaks will be equal to the total number of unique stereoisomers in the mixture.
    *   Total unique stereoisomers = 2 (from Reaction I) + 3 (from Reaction II) = **5 stereoisomers**.
    *   A chiral column is capable of resolving all five of these distinct stereoisomers.
    *   **Total Peaks:** **5 peaks**.

### **Final Conclusion**

The analysis shows that the chemist will observe **5 peaks in the chiral HPLC** and **4 peaks in the normal-phase HPLC**. This result corresponds to option B.

<<<B>>>
"""

# Run the check
result = check_correctness_of_chemistry_answer(question_full_text, final_answer_from_llm)
print(result)