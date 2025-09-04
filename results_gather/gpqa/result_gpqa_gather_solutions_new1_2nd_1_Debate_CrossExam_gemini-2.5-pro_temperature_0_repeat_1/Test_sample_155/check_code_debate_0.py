import re

def check_answer_correctness(final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the chemistry HPLC problem.

    Args:
        final_answer_text: The text of the final answer, including the <<<...>>> tag.

    Returns:
        "Correct" if the answer is correct, or a string explaining the error.
    """

    # --- Step 1: Define Chemical and Chromatographic Principles ---

    # Stereochemical rules for the anti-dihydroxylation of a symmetrical alkene.
    # The product types are 'meso' or 'enantiomer_pair'.
    def get_reaction_products(alkene_type: str) -> list:
        if alkene_type == 'trans':  # (E)-alkene + anti-addition -> meso
            return ['meso']
        elif alkene_type == 'cis':  # (Z)-alkene + anti-addition -> racemic (enantiomeric pair)
            return ['enantiomer_pair']
        return []

    # HPLC separation rules.
    def count_standard_hplc_peaks(product_types: list) -> int:
        """Calculates peaks on an achiral column where enantiomers co-elute."""
        # Use a set to count unique peak types.
        # The meso compound is a diastereomer to the enantiomeric pair, so they separate.
        # The two enantiomers within the pair do not separate from each other.
        peak_types = set(product_types)
        return len(peak_types)

    def count_chiral_hplc_peaks(product_types: list) -> int:
        """Calculates peaks on a chiral column where all stereoisomers are resolved."""
        # A meso compound gives 1 peak.
        # An enantiomeric pair is resolved into 2 peaks.
        peaks = 0
        if 'meso' in product_types:
            peaks += 1
        if 'enantiomer_pair' in product_types:
            peaks += 2
        return peaks

    # --- Step 2: Apply Principles to the Problem ---

    # Reaction 1: (E)-oct-4-ene is a 'trans' alkene.
    products_r1 = get_reaction_products('trans')

    # Reaction 2: (Z)-oct-4-ene is a 'cis' alkene.
    products_r2 = get_reaction_products('cis')

    # The final mixture contains the products of both reactions.
    combined_product_types = products_r1 + products_r2
    # This results in ['meso', 'enantiomer_pair']

    # Calculate the theoretically correct number of peaks.
    expected_standard_peaks = count_standard_hplc_peaks(combined_product_types)
    expected_chiral_peaks = count_chiral_hplc_peaks(combined_product_types)

    # --- Step 3: Parse and Evaluate the LLM's Final Answer ---

    # Define the multiple-choice options from the question.
    options = {
        'A': {'standard': 4, 'chiral': 4},
        'B': {'standard': 3, 'chiral': 4},
        'C': {'standard': 2, 'chiral': 2},
        'D': {'standard': 2, 'chiral': 3}
    }

    # Extract the chosen option letter from the final answer text.
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not find a valid answer tag like '<<<A>>>' in the provided text."

    chosen_option_letter = match.group(1)
    chosen_option_values = options[chosen_option_letter]

    # --- Step 4: Compare Expected vs. Chosen Answer ---

    if (chosen_option_values['standard'] == expected_standard_peaks and
            chosen_option_values['chiral'] == expected_chiral_peaks):
        return "Correct"
    else:
        # Provide a detailed reason for the error.
        reason = (
            f"Incorrect. The final answer chose option '{chosen_option_letter}', which corresponds to "
            f"{chosen_option_values['standard']} standard peaks and {chosen_option_values['chiral']} chiral peaks.\n\n"
            f"The correct analysis is:\n"
            f"1. The final mixture contains three stereoisomers: one meso compound and one pair of enantiomers.\n"
            f"2. On a standard (achiral) HPLC, diastereomers separate but enantiomers co-elute. This results in 2 peaks (one for the meso compound, one for the enantiomeric pair).\n"
            f"3. On a chiral HPLC, all stereoisomers are resolved. This results in 3 peaks (one for the meso compound, and two for the separated enantiomers).\n\n"
            f"Therefore, the correct answer is {expected_standard_peaks} peaks in standard HPLC and {expected_chiral_peaks} peaks in chiral HPLC, which corresponds to option D."
        )
        return reason

# The user provided a final answer block to check.
final_answer_from_user = """
Here is a step-by-step analysis of the problem, synthesizing the reasoning from the provided answers and identifying the key points required for the correct solution.

### **Step 1: Analyze the Chemical Reaction**

The core of the problem is to identify the products of the two reactions. The described two-step process, treatment with mCPBA followed by aqueous acid (H₃O⁺), results in the **anti-dihydroxylation** of an alkene. This means two hydroxyl (-OH) groups are added across the double bond, one to each face.

### **Step 2: Determine the Products of Reaction 1**

*   **Starting Material:** (E)-oct-4-ene, which is a symmetrical *trans*-alkene.
*   **Transformation:** Anti-dihydroxylation.
*   **Stereochemical Rule:** The anti-addition of two identical groups to a symmetrical *trans*-alkene results in a **meso compound**. A meso compound is achiral despite having chiral centers due to an internal plane of symmetry.
*   **Product of Reaction 1:** A single product is formed: **meso-octane-4,5-diol**.

### **Step 3: Determine the Products of Reaction 2**

*   **Starting Material:** (Z)-oct-4-ene, which is a symmetrical *cis*-alkene.
*   **Transformation:** Anti-dihydroxylation.
*   **Stereochemical Rule:** The anti-addition of two identical groups to a symmetrical *cis*-alkene results in a **racemic mixture** (a 1:1 mixture of enantiomers).
*   **Product of Reaction 2:** A pair of enantiomers is formed: **(4R,5R)-octane-4,5-diol** and **(4S,5S)-octane-4,5-diol**.

### **Step 4: Analyze the Combined Product Mixture**

The chemist combines the products from both reactions. The final mixture contains three distinct stereoisomers:
1.  **meso-octane-4,5-diol** (from Reaction 1)
2.  **(4R,5R)-octane-4,5-diol** (from Reaction 2)
3.  **(4S,5S)-octane-4,5-diol** (from Reaction 2)

The relationships between these isomers are critical:
*   The (4R,5R) and (4S,5S) compounds are **enantiomers** of each other.
*   The meso compound is a **diastereomer** of both the (4R,5R) and (4S,5S) compounds.

### **Step 5: Predict the Chromatographic Results**

**1. Standard (achiral) Reverse-Phase HPLC:**
*   **Principle:** This technique separates compounds based on differences in physical properties, such as polarity. It cannot distinguish between enantiomers.
*   **Analysis:**
    *   Enantiomers have identical physical properties in an achiral environment. Therefore, the (4R,5R)-diol and (4S,5S)-diol will not be separated and will elute together as a single peak.
    *   Diastereomers have different physical properties. The meso-diol is a diastereomer of the enantiomeric pair and will have a different retention time, eluting as a separate peak.
*   **Conclusion:** The chromatogram will show **2 peaks**.

**2. Chiral HPLC:**
*   **Principle:** This technique uses a chiral stationary phase, which creates a chiral environment. This allows the column to interact differently with each enantiomer, enabling their separation. It also separates diastereomers.
*   **Analysis:** The chiral column can distinguish between all three unique stereoisomers in the mixture.
    *   The meso-diol will give one peak.
    *   The (4R,5R)-enantiomer will give a second peak.
    *   The (4S,5S)-enantiomer will give a third peak.
*   **Conclusion:** The chromatogram will show **3 peaks**.

### **Careful Points to Consider**

*   **Reaction Stereochemistry:** It is essential to correctly identify the overall reaction as an *anti*-dihydroxylation and apply the correct stereochemical rules: *trans*-alkene + *anti*-addition → *meso* compound; *cis*-alkene + *anti*-addition → *racemic* mixture.
*   **Isomer Relationships:** Correctly identifying the products as one meso compound and one pair of enantiomers (which are diastereomers of the meso compound) is crucial.
*   **HPLC Principles:** One must clearly distinguish between the capabilities of achiral and chiral HPLC. Achiral HPLC separates diastereomers but not enantiomers. Chiral HPLC separates both.
*   **Matching Analysis to Options:** After determining the outcome (2 peaks in standard HPLC and 3 peaks in chiral HPLC), it is essential to select the correct corresponding option from the list.

The analysis consistently shows there will be 2 peaks in the standard HPLC and 3 peaks in the chiral HPLC. This corresponds to option D.

<<<D>>>
"""

# Run the check
result = check_answer_correctness(final_answer_from_user)
print(result)