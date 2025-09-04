import re

def check_correctness(question, final_answer_text):
    """
    Checks the correctness of the final answer based on chemical principles.

    Args:
        question (str): The text of the question, including the options.
        final_answer_text (str): The text of the final answer, including the <<<B>>> part.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # Step 1: Determine the correct number of peaks based on chemical principles.
    
    # Reaction I: (S)-5-methoxyhexan-3-one reduction
    # A chiral starting material with one stereocenter (C5) is reduced at a ketone (C3),
    # creating a new stereocenter. This results in two products: (3R, 5S) and (3S, 5S).
    # These products are diastereomers.
    # Number of products from Reaction I = 2 (diastereomers)
    rxn1_products = {"type": "diastereomers", "count": 2}

    # Reaction II: Pentane-2,4-dione reduction
    # An achiral, symmetric starting material is reduced at two ketones, creating two new
    # stereocenters (C2, C4). This results in a pair of enantiomers ((2R,4R) and (2S,4S))
    # and a meso compound ((2R,4S)).
    # Number of unique stereoisomers from Reaction II = 3
    rxn2_products = {"type": "enantiomers_and_meso", "enantiomer_pairs": 1, "meso_compounds": 1}

    # Calculate peaks for Normal-Phase (achiral) HPLC
    # It separates diastereomers but not enantiomers.
    # Products from Rxn I and Rxn II are structurally different and will not co-elute.
    # Rxn I: 2 diastereomers -> 2 peaks
    # Rxn II: 1 enantiomeric pair (1 peak) + 1 meso compound (1 peak) -> 2 peaks
    calculated_normal_peaks = rxn1_products["count"] + rxn2_products["enantiomer_pairs"] + rxn2_products["meso_compounds"]

    # Calculate peaks for Chiral HPLC
    # It separates all unique stereoisomers.
    # Rxn I: 2 stereoisomers -> 2 peaks
    # Rxn II: 3 stereoisomers (2 enantiomers + 1 meso) -> 3 peaks
    calculated_chiral_peaks = rxn1_products["count"] + (rxn2_products["enantiomer_pairs"] * 2) + rxn2_products["meso_compounds"]

    # Step 2: Parse the options and the selected answer from the provided text.
    
    def parse_options(text):
        options = {}
        option_lines = re.findall(r"^[A-D]\) .*", text, re.MULTILINE)
        for line in option_lines:
            key = line[0]
            match_both = re.search(r"(\d+) peaks in both", line)
            if match_both:
                peaks = int(match_both.group(1))
                options[key] = {"chiral": peaks, "normal": peaks}
                continue
            
            match_separate = re.search(r"(\d+)\s+Peaks in chiral HPLC and (\d+)\s+peaks in normal-phase HPLC", line, re.IGNORECASE)
            if match_separate:
                chiral_peaks = int(match_separate.group(1))
                normal_peaks = int(match_separate.group(2))
                options[key] = {"chiral": chiral_peaks, "normal": normal_peaks}
        return options

    options_map = parse_options(question)
    if not options_map:
        return "Error: Could not parse the options from the question text."

    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not find the final answer in the format <<<X>>>."
    
    selected_option_key = match.group(1)
    
    if selected_option_key not in options_map:
        return f"Error: The selected option '{selected_option_key}' is not a valid choice."

    claimed_peaks = options_map[selected_option_key]

    # Step 3: Compare the calculated results with the claimed answer.
    
    if calculated_chiral_peaks != claimed_peaks["chiral"]:
        return (f"Incorrect. The number of peaks in chiral HPLC is wrong. "
                f"Calculation: Reaction I produces 2 diastereomers. Reaction II produces 3 unique stereoisomers (an enantiomeric pair and a meso compound). "
                f"A chiral column separates all 5 unique stereoisomers. "
                f"Expected chiral peaks: {calculated_chiral_peaks}. The answer claims: {claimed_peaks['chiral']}.")

    if calculated_normal_peaks != claimed_peaks["normal"]:
        return (f"Incorrect. The number of peaks in normal-phase HPLC is wrong. "
                f"Calculation: In an achiral column, the 2 diastereomers from Reaction I give 2 peaks. The enantiomeric pair from Reaction II gives 1 peak, and the meso compound gives 1 peak. "
                f"Total expected normal-phase peaks: 2 + 1 + 1 = {calculated_normal_peaks}. The answer claims: {claimed_peaks['normal']}.")

    return "Correct"

# The question and final answer text are provided in the context.
question_text = """
An organic chemist performs two reactions.
Reaction I:
(S)-5-methoxyhexan-3-one is treated with LAH followed by an acidic workup.
Reaction II:
Pentane-2,4-dione is treated with excess NaBH4 followed by an acidic workup.
The progress of the reactions is monitored using TLC. After the reactions are 100% completed, the products from both reactions are combined, and run on both a normal-phase HPLC column and an HPLC column packed with the chiral stationary phase.
How many peaks does the chemist observe in the chromatograms?

A) 4 Peaks in chiral HPLC and 2 peaks in normal-phase HPLC
B) 5 Peaks in chiral HPLC and 4 peaks in normal-phase HPLC
C) 3 peaks in both chiral HPLC and normal-phase HPLC
D) 3 Peaks in chiral HPLC and 2 peaks in normal-phase HPLC
"""

final_answer_text = """
Here is a step-by-step analysis to determine the final answer, taking into account the reasoning from the provided candidate answers.

### **Step 1: Analyze the Products of Reaction I**

*   **Starting Material:** (S)-5-methoxyhexan-3-one. This molecule is chiral and has a single pre-existing stereocenter at the C5 position with an (S) configuration.
*   **Reaction:** The ketone group at the C3 position is reduced to a secondary alcohol using Lithium Aluminum Hydride (LAH). This reduction creates a new stereocenter at C3.
*   **Products:** The original stereocenter at C5 is not involved in the reaction and retains its (S) configuration. The reduction of the planar carbonyl group can occur from two faces, leading to both (R) and (S) configurations at the new C3 stereocenter. This results in two products:
    1.  (3R, 5S)-5-methoxyhexan-3-ol
    2.  (3S, 5S)-5-methoxyhexan-3-ol
*   **Relationship:** These two products are **diastereomers** because they have two stereocenters, one of which is the same (5S) and one of which is different (3R vs. 3S).
*   **Conclusion for Reaction I:** Two distinct diastereomeric products are formed.

### **Step 2: Analyze the Products of Reaction II**

*   **Starting Material:** Pentane-2,4-dione. This is a symmetric and achiral molecule.
*   **Reaction:** Both ketone groups are reduced to secondary alcohols using excess Sodium Borohydride (NaBH₄). This creates two new stereocenters at the C2 and C4 positions.
*   **Products:** Since the starting material and reagent are achiral, all possible stereoisomers of the product (pentane-2,4-diol) are formed.
    1.  (2R, 4R)-pentane-2,4-diol
    2.  (2S, 4S)-pentane-2,4-diol
    3.  (2R, 4S)-pentane-2,4-diol
*   **Relationship:**
    *   The (2R, 4R) and (2S, 4S) isomers are non-superimposable mirror images, making them a pair of **enantiomers**. They are produced as a 1:1 racemic mixture.
    *   The (2R, 4S) isomer has an internal plane of symmetry, which makes it an achiral **meso compound**.
    *   The meso compound is a diastereomer of the enantiomeric pair.
*   **Conclusion for Reaction II:** Three distinct stereoisomers are formed.

### **Step 3: Analyze the Combined Mixture with Normal-Phase HPLC**

*   **Principle:** A normal-phase HPLC column is achiral. It separates compounds based on differences in physical properties like polarity. It can separate constitutional isomers and diastereomers, but it **cannot** separate enantiomers, as they have identical physical properties in an achiral environment.
*   **Peak Count:**
    *   The two diastereomers from Reaction I have different physical properties and will be separated. This gives **2 peaks**.
    *   From Reaction II, the enantiomeric pair ((2R, 4R) and (2S, 4S)) will elute together as **1 peak**.
    *   The meso compound from Reaction II is a diastereomer of the enantiomeric pair, has different physical properties, and will elute as a separate, **second peak**.
    *   The products from Reaction I (C7 ether-alcohols) and Reaction II (C5 diols) are structurally different and will not accidentally elute at the same time.
*   **Total Peaks in Normal-Phase HPLC:** 2 (from Reaction I) + 2 (from Reaction II) = **4 peaks**.

### **Step 4: Analyze the Combined Mixture with Chiral HPLC**

*   **Principle:** A chiral HPLC column has a chiral stationary phase that interacts differently with all stereoisomers, allowing for the separation of **both diastereomers and enantiomers**.
*   **Peak Count:**
    *   The column is capable of resolving every unique stereoisomer in the mixture.
    *   Total unique stereoisomers = 2 (from Reaction I) + 3 (from Reaction II) = **5 stereoisomers**.
*   **Total Peaks in Chiral HPLC:** **5 peaks**.

### **Final Conclusion**

The detailed analysis shows that the chemist will observe **5 peaks in the chiral HPLC** and **4 peaks in the normal-phase HPLC**.

Now, we match this result to the given options:
A) 4 Peaks in chiral HPLC and 2 peaks in normal-phase HPLC
B) 5 Peaks in chiral HPLC and 4 peaks in normal-phase HPLC
C) 3 peaks in both chiral HPLC and normal-phase HPLC
D) 3 Peaks in chiral HPLC and 2 peaks in normal-phase HPLC

The correct option is B. The vast majority of the candidate answers reached the same conclusion regarding the number of peaks, even if some incorrectly mapped their reasoning to the final letter choice.

<<<B>>>
"""

# Run the check
result = check_correctness(question_text, final_answer_text)
print(result)