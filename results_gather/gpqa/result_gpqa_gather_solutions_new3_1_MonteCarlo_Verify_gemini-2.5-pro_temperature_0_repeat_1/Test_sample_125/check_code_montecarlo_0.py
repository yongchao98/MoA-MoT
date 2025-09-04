import re

def check_correctness(question: str, llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer to the organic chemistry HPLC question.

    The function will codify the chemical logic:
    1. Determine the number and type of products from Reaction I.
    2. Determine the number and type of products from Reaction II.
    3. Calculate the expected number of peaks for each HPLC type based on stereochemical principles.
    4. Compare the calculated result with the provided answer.
    """

    # --- Step 1: Analyze Reaction I ---
    # Reactant: (S)-5-methoxyhexan-3-one (chiral, 1 stereocenter)
    # Reaction: Reduction of ketone at C3 creates a new stereocenter.
    # Products: (3R, 5S)- and (3S, 5S)-5-methoxyhexan-3-ol.
    # Relationship: These are diastereomers.
    # Normal HPLC: Diastereomers have different physical properties and will separate. -> 2 peaks.
    # Chiral HPLC: Diastereomers are different compounds and will separate. -> 2 peaks.
    peaks_from_rxn1_normal = 2
    peaks_from_rxn1_chiral = 2

    # --- Step 2: Analyze Reaction II ---
    # Reactant: Pentane-2,4-dione (achiral)
    # Reaction: Reduction of both ketones creates two new stereocenters (C2, C4).
    # Products: (2R,4R)-, (2S,4S)-, and meso-(2R,4S)-pentane-2,4-diol.
    # Relationship: One enantiomeric pair and one meso compound. The meso compound is a diastereomer of the enantiomeric pair.
    # Normal HPLC: The enantiomeric pair has identical properties in an achiral environment and co-elutes (1 peak). The meso compound is a diastereomer and separates (1 peak). -> 2 peaks total.
    # Chiral HPLC: The enantiomeric pair is resolved by the chiral column (2 peaks). The meso compound is a distinct stereoisomer and also gives a peak (1 peak). -> 3 peaks total.
    peaks_from_rxn2_normal = 2
    peaks_from_rxn2_chiral = 3

    # --- Step 3: Calculate total expected peaks ---
    # We assume that the products from Reaction I (C7 ether-alcohols) and Reaction II (C5 diols)
    # are structurally different enough to be resolved from each other in both systems.
    expected_normal_peaks = peaks_from_rxn1_normal + peaks_from_rxn2_normal
    expected_chiral_peaks = peaks_from_rxn1_chiral + peaks_from_rxn2_chiral

    # --- Step 4: Parse the LLM's final answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The final answer is not in the required format '<<<A>>>', '<<<B>>>', etc."

    llm_choice = match.group(1)

    # --- Step 5: Define the options from the question ---
    options = {
        'A': {'chiral': 5, 'normal': 4},
        'B': {'chiral': 3, 'normal': 2},
        'C': {'chiral': 3, 'normal': 3},
        'D': {'chiral': 4, 'normal': 2}
    }

    if llm_choice not in options:
        return f"Incorrect: The chosen option '{llm_choice}' is not a valid option (A, B, C, or D)."

    llm_peaks = options[llm_choice]
    llm_chiral_peaks = llm_peaks['chiral']
    llm_normal_peaks = llm_peaks['normal']

    # --- Step 6: Compare expected results with the LLM's answer ---
    if llm_chiral_peaks == expected_chiral_peaks and llm_normal_peaks == expected_normal_peaks:
        return "Correct"
    else:
        correct_choice = 'A' # Based on our calculation
        reason = (f"Incorrect. The chosen answer is {llm_choice}, but the correct answer is {correct_choice}.\n"
                  f"Reasoning:\n"
                  f"1.  **Reaction I** produces 2 diastereomers. These give 2 peaks in both normal and chiral HPLC.\n"
                  f"2.  **Reaction II** produces 3 stereoisomers (an enantiomeric pair and a meso compound). These give 2 peaks in normal HPLC (enantiomers co-elute) and 3 peaks in chiral HPLC (enantiomers separate).\n"
                  f"3.  **Total Normal Peaks** = {peaks_from_rxn1_normal} (from Rxn I) + {peaks_from_rxn2_normal} (from Rxn II) = {expected_normal_peaks}.\n"
                  f"4.  **Total Chiral Peaks** = {peaks_from_rxn1_chiral} (from Rxn I) + {peaks_from_rxn2_chiral} (from Rxn II) = {expected_chiral_peaks}.\n"
                  f"Therefore, the correct result is {expected_chiral_peaks} peaks in chiral HPLC and {expected_normal_peaks} peaks in normal-phase HPLC, which corresponds to option {correct_choice}.")
        return reason

# The user provided the question and a final answer block to check.
question_text = """
An organic chemist performs two reactions.
Reaction I:
(S)-5-methoxyhexan-3-one is treated with LAH followed by an acidic workup.
Reaction II:
Pentane-2,4-dione is treated with excess NaBH4 followed by an acidic workup.
The progress of the reactions is monitored using TLC. After the reactions are 100% completed, the products from both reactions are combined, and run on both a normal-phase HPLC column and an HPLC column packed with the chiral stationary phase.
How many peaks does the chemist observe in the chromatograms?

A) 5 Peaks in chiral HPLC and 4 peaks in normal-phase HPLC
B) 3 Peaks in chiral HPLC and 2 peaks in normal-phase HPLC
C) 3 peaks in both chiral HPLC and normal-phase HPLC
D) 4 Peaks in chiral HPLC and 2 peaks in normal-phase HPLC
"""

llm_final_answer = """
Here is a step-by-step analysis to determine the correct answer.

### **Careful Points to Consider**

1.  **Stereochemistry of Starting Materials:** Note that the reactant in Reaction I is chiral, while the reactant in Reaction II is achiral.
2.  **Creation of New Stereocenters:** Both reactions involve the reduction of ketones, which creates new stereocenters.
3.  **Stereochemical Outcome:** The relationship between the products (enantiomers, diastereomers, meso compounds) is crucial.
4.  **Principles of HPLC Separation:**
    *   **Normal-Phase HPLC (Achiral Column):** Separates compounds based on differences in physical properties like polarity. It can separate diastereomers but **cannot** separate enantiomers.
    *   **Chiral HPLC (Chiral Column):** Uses a chiral stationary phase to interact differently with enantiomers, allowing for their separation. It can separate all unique stereoisomers.
5.  **Combined Mixture:** The final analysis is performed on the combined products of both reactions. We assume that the products from the two different reactions are structurally distinct enough to be separated from each other.

---

### **Step 1: Analyze the Products of Reaction I**

*   **Starting Material:** (S)-5-methoxyhexan-3-one. This is a chiral molecule with one existing stereocenter at the C5 position.
*   **Reaction:** The ketone at C3 is reduced to a secondary alcohol by LAH. This reduction creates a new stereocenter at C3.
*   **Products:** The original stereocenter at C5 is unaffected. The new stereocenter at C3 can be either (R) or (S). This results in the formation of two products:
    1.  (3R, 5S)-5-methoxyhexan-3-ol
    2.  (3S, 5S)-5-methoxyhexan-3-ol
*   **Relationship:** These two products are **diastereomers** because one stereocenter is the same (5S) while the other is different (3R vs. 3S).

**Conclusion for Reaction I:** Two diastereomeric products are formed.

### **Step 2: Analyze the Products of Reaction II**

*   **Starting Material:** Pentane-2,4-dione. This is a symmetric, achiral molecule.
*   **Reaction:** Both ketone groups are reduced to secondary alcohols by excess NaBH₄, creating two new stereocenters at C2 and C4.
*   **Products:** The reaction produces all possible stereoisomers of pentane-2,4-diol.
    1.  (2R, 4R)-pentane-2,4-diol
    2.  (2S, 4S)-pentane-2,4-diol
    3.  (2R, 4S)-pentane-2,4-diol (which is identical to 2S, 4R)
*   **Relationship:**
    *   The (2R, 4R) and (2S, 4S) isomers are a pair of **enantiomers**.
    *   The (2R, 4S) isomer has an internal plane of symmetry, making it an achiral **meso compound**.
    *   The meso compound is a diastereomer of the enantiomeric pair.

**Conclusion for Reaction II:** Three distinct stereoisomers are formed (a pair of enantiomers and a meso compound).

### **Step 3: Analyze the Combined Mixture with HPLC**

The final mixture contains a total of five distinct stereoisomers (two from Reaction I and three from Reaction II).

#### **Normal-Phase HPLC (Achiral Column)**

This column separates diastereomers but not enantiomers.
*   **From Reaction I:** The two diastereomers will be separated. This gives **2 peaks**.
*   **From Reaction II:** The enantiomeric pair will elute together as one peak. The meso compound, being a diastereomer of the pair, will elute separately. This gives **2 peaks**.
*   **Total Peaks:** The products from Reaction I (C7 ether-alcohols) and Reaction II (C5 diols) are structurally different and will be resolved. Total peaks = 2 (from Rxn I) + 2 (from Rxn II) = **4 peaks**.

#### **Chiral HPLC**

This column separates all unique stereoisomers, including enantiomers.
*   **From Reaction I:** The two diastereomers will be separated. This gives **2 peaks**.
*   **From Reaction II:** The two enantiomers will be resolved into two separate peaks. The meso compound will also give its own distinct peak. This gives **3 peaks**.
*   **Total Peaks:** Total peaks = 2 (from Rxn I) + 3 (from Rxn II) = **5 peaks**.

### **Final Conclusion**

The analysis predicts **5 peaks in chiral HPLC** and **4 peaks in normal-phase HPLC**. This corresponds to option A.

<<<A>>>
"""

# Run the check
result = check_correctness(question_text, llm_final_answer)
print(result)