import re

def check_chemistry_answer(llm_answer):
    """
    Checks the correctness of the answer for the given chemistry problem.

    The function verifies the chosen chemical formula against constraints derived
    from spectroscopic data (FTIR, 1H NMR, Mass Spec).

    Constraints:
    1. Must contain 2 oxygen atoms (due to carboxylic acid).
    2. Must have a Degree of Unsaturation (DoU) of 2 (for C=O and C=C).
    """

    # Extract the letter from the answer format, e.g., 'A' from '<<<A>>>'
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "Incorrect. The answer format is invalid. It should be like '<<<A>>>'."
    
    chosen_option = match.group(1)

    # Define the options from the question
    options = {
        'A': 'C6H10O2',
        'B': 'C6H10O',
        'C': 'C6H12O2',
        'D': 'C6H12O'
    }
    
    # Get the formula corresponding to the chosen answer
    chosen_formula = options.get(chosen_option)
    if not chosen_formula:
        # This case should not be hit if the regex is correct, but it's good practice
        return f"Invalid option '{chosen_option}' selected."

    # Helper function to parse a chemical formula string
    def parse_formula(formula_str):
        """Parses a formula like 'C6H10O2' into a dictionary of element counts."""
        counts = {}
        # Find all element-number pairs (e.g., C6, H10, O2)
        element_matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula_str)
        for element, count in element_matches:
            # If count is empty, it means 1 atom (e.g., in H2O)
            counts[element] = int(count) if count else 1
        return counts

    # Parse the chosen formula
    atom_counts = parse_formula(chosen_formula)
    
    # --- Constraint 1: Check for 2 Oxygen atoms ---
    # Evidence: FTIR (broad 3000, strong 1700) and Mass Spec (m/z=45) indicate a carboxylic acid.
    required_oxygens = 2
    actual_oxygens = atom_counts.get('O', 0)
    if actual_oxygens != required_oxygens:
        return (f"Incorrect. The answer '{chosen_option}' corresponds to {chosen_formula}. "
                f"The evidence points to a carboxylic acid, which requires {required_oxygens} oxygen atoms, "
                f"but this formula has {actual_oxygens}.")

    # --- Constraint 2: Check Degree of Unsaturation (DoU) ---
    # Evidence: Carboxylic acid (C=O) and alkene (C=C) require a DoU of 2.
    required_dou = 2
    num_carbons = atom_counts.get('C', 0)
    num_hydrogens = atom_counts.get('H', 0)
    
    # DoU formula for CxHyOz: C + 1 - (H/2)
    actual_dou = num_carbons + 1 - (num_hydrogens / 2)
    
    if actual_dou != required_dou:
        return (f"Incorrect. The answer '{chosen_option}' corresponds to {chosen_formula}. "
                f"The evidence points to one C=O and one C=C bond, requiring a Degree of Unsaturation (DoU) of {required_dou}. "
                f"The calculated DoU for this formula is {actual_dou}.")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# The final answer provided by the LLM
final_answer_from_llm = """
Based on a step-by-step analysis of the spectroscopic data and the candidate answers, here is the definitive resolution.

### Step-by-Step Analysis

The core of this problem is to interpret the provided spectroscopic data to establish a set of rules or constraints that the correct molecular formula must obey.

**1. Analysis of the FTIR Spectrum**
*   **Very broad peak at 3000 cm⁻¹:** This is a classic and unambiguous signature of the O-H stretching vibration in a **carboxylic acid (-COOH)**. The extreme broadness is caused by strong hydrogen bonding.
*   **Strong peak at 1700 cm⁻¹:** This peak is in the carbonyl (C=O) stretching region. Its presence, combined with the broad O-H peak, confirms the carboxylic acid functional group.
*   **Strong peak at 1650 cm⁻¹:** This peak is characteristic of a carbon-carbon double bond (C=C) stretch, indicating the molecule is an **alkene**.

**2. Analysis of the ¹H NMR Spectrum**
*   **Peaks for vinyl-hydrogens:** Vinyl hydrogens are protons directly attached to the carbons of a C=C double bond. This observation provides direct confirmation of the alkene functional group, supporting the FTIR data.

**3. Analysis of the Mass Spectrometry Data**
*   **Fragment peak at m/z = 45:** This is a highly diagnostic fragment for a carboxylic acid. It corresponds to the mass of the carboxyl fragment, [COOH]⁺. (Mass = 12 for C + 16 for O + 16 for O + 1 for H = 45).

### Careful Points and Constraints

The combined evidence from all three techniques points to a molecule with two specific functional groups. This allows us to establish firm constraints.

*   **Careful Point 1: Oxygen Count.** The evidence for a carboxylic acid (-COOH) is overwhelming. This functional group contains **two oxygen atoms**. This constraint immediately eliminates any formulas with only one oxygen atom.
    *   Eliminated: B) C6H10O and D) C6H12O.
    *   Remaining Candidates: A) C6H10O2 and C) C6H12O2.

*   **Careful Point 2: Degree of Unsaturation (DoU).** The Degree of Unsaturation represents the total number of rings and pi (π) bonds in a molecule. The identified functional groups require a specific DoU:
    *   The C=O double bond in the carboxylic acid accounts for **one** degree of unsaturation.
    *   The C=C double bond in the alkene accounts for a **second** degree of unsaturation.
    *   Therefore, the correct molecule must have a total DoU of **2**.

*   **Careful Point 3: Evaluating Remaining Candidates.** We now calculate the DoU for the two remaining formulas using the formula `DoU = C + 1 - (H/2)`.
    *   **For A) C6H10O2:** DoU = 6 + 1 - (10/2) = 7 - 5 = **2**. This formula perfectly matches all constraints.
    *   **For C) C6H12O2:** DoU = 6 + 1 - (12/2) = 7 - 6 = **1**. This formula is incorrect. Its DoU of 1 can only account for the carboxylic acid's C=O bond, but not the additional C=C bond confirmed by FTIR and NMR.

### Final Conclusion

The chemical reasoning is unambiguous and consistently supported by all provided candidate answers: the only formula that fits all the spectroscopic evidence is **C6H10O2**.

The final step is to map this correct formula to the options given in the question header:
A) C6H10O2
B) C6H10O
C) C6H12O2
D) C6H12O

The correct formula, C6H10O2, corresponds to option A. Many of the candidate answers performed the chemical analysis correctly but made a clerical error in mapping the formula to the final letter. The correct reasoning leads to A.

<<<A>>>
"""

# Run the check
result = check_chemistry_answer(final_answer_from_llm)
print(result)