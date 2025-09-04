import re

def check_answer(final_answer_text):
    """
    Checks if the final answer correctly identifies the chemical formula based on spectroscopic data.

    The function verifies two main things:
    1. The chosen formula must have 2 oxygen atoms (for the carboxylic acid).
    2. The chosen formula must have a Degree of Unsaturation (DoU) of 2 (for C=O and C=C bonds).
    """

    # --- Step 1: Define constraints from the problem description ---
    required_oxygens = 2
    required_dou = 2

    # --- Step 2: Define helper functions ---
    def parse_formula(formula_str):
        """Parses a chemical formula string into a dictionary of atom counts."""
        try:
            c = int(re.search(r'C(\d+)', formula_str).group(1))
            h = int(re.search(r'H(\d+)', formula_str).group(1))
            o_match = re.search(r'O(\d*)', formula_str)
            o = int(o_match.group(1)) if o_match and o_match.group(1) else (1 if o_match else 0)
            return {'C': c, 'H': h, 'O': o}
        except (AttributeError, ValueError):
            return None

    def calculate_dou(atoms):
        """Calculates the Degree of Unsaturation for a given set of atoms."""
        # Formula: DoU = C + 1 - (H/2)
        return atoms['C'] + 1 - (atoms['H'] / 2)

    # --- Step 3: Extract the options and the selected answer from the text ---
    # The final answer being checked provides its own mapping of letters to formulas.
    # We must use this specific mapping to check its internal consistency.
    options_map = {}
    # Example line to match: "A) C6H10O" or "*   **A) C6H10O2:**"
    option_lines = re.findall(r'([A-D])\)\s+(C\d+H\d+O\d*)', final_answer_text)
    if not option_lines:
        return "Incorrect. The answer text does not clearly list the options A, B, C, D with their formulas."
        
    for letter, formula in option_lines:
        options_map[letter] = formula

    # Extract the final choice, e.g., 'C' from '<<<C>>>'
    final_choice_match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not final_choice_match:
        return "Incorrect. The final answer is not provided in the required '<<<X>>>' format."
    
    selected_letter = final_choice_match.group(1)
    
    if selected_letter not in options_map:
        return f"Incorrect. The selected answer '{selected_letter}' does not correspond to any option listed in the reasoning."

    selected_formula = options_map[selected_letter]

    # --- Step 4: Verify the selected formula against the constraints ---
    atoms = parse_formula(selected_formula)
    if not atoms:
        return f"Incorrect. The selected formula '{selected_formula}' could not be parsed."

    # Check oxygen count
    if atoms['O'] != required_oxygens:
        return (f"Incorrect. The selected formula {selected_formula} has {atoms['O']} oxygen atoms, "
                f"but the presence of a carboxylic acid requires {required_oxygens} oxygen atoms.")

    # Check Degree of Unsaturation
    dou = calculate_dou(atoms)
    if dou != required_dou:
        return (f"Incorrect. The selected formula {selected_formula} has a Degree of Unsaturation (DoU) of {int(dou)}, "
                f"but the presence of both a C=O and a C=C bond requires a DoU of {required_dou}.")

    # --- Step 5: Verify that the chosen answer is the *only* correct one ---
    correct_options = []
    for letter, formula in options_map.items():
        f_atoms = parse_formula(formula)
        if f_atoms and f_atoms['O'] == required_oxygens and calculate_dou(f_atoms) == required_dou:
            correct_options.append(letter)
    
    if len(correct_options) > 1:
        return f"Incorrect. The question is ambiguous as options {correct_options} both satisfy the criteria."
    
    if selected_letter not in correct_options:
        return (f"Incorrect. The selected option {selected_letter} ({selected_formula}) does not meet the criteria. "
                f"The correct option is {correct_options[0]} ({options_map[correct_options[0]]}).")

    return "Correct"

# The final answer provided by the user
final_answer_text = """
Here is a step-by-step analysis to determine the final answer.

**Plan:**
1.  Break down the spectroscopic data (FTIR, NMR, MS) to identify the key functional groups present in the unknown compound.
2.  Use the identified functional groups to establish mandatory criteria for the molecular formula (e.g., number of oxygen atoms, degree of unsaturation).
3.  Evaluate each of the four candidate formulas against these criteria to find the only one that matches.
4.  Carefully map the correct formula to its corresponding letter (A, B, C, or D) to select the final answer.

**Execution:**

**Step 1: Analyze the Spectroscopic Data**

*   **FTIR Spectrum:**
    *   A **very broad peak at 3000 cm⁻¹** combined with a **strong peak at 1700 cm⁻¹** is the definitive signature of a **carboxylic acid (-COOH)** functional group. The broadness is due to O-H hydrogen bonding, and 1700 cm⁻¹ is the C=O (carbonyl) stretch.
    *   A **strong peak at 1650 cm⁻¹** is characteristic of a **carbon-carbon double bond (C=C)**, indicating the molecule is an alkene.
*   **¹H NMR Spectrum:**
    *   The presence of **vinyl-hydrogens** explicitly confirms the existence of the C=C double bond, supporting the FTIR data.
*   **Mass Spectrometry:**
    *   A **fragment peak at m/z = 45** is a classic and highly diagnostic fragment for a carboxylic acid. It corresponds to the mass of the [COOH]⁺ ion (12 for C + 16 for O + 16 for O + 1 for H = 45).

**Conclusion from Data:** The combined evidence from all three techniques proves that the unknown compound is an **unsaturated carboxylic acid**.

**Step 2: Establish Criteria for the Molecular Formula**

Based on the conclusion from Step 1, the correct molecular formula must satisfy two critical conditions:
1.  **Oxygen Count:** It must contain a carboxylic acid group (-COOH), which means the formula must have exactly **two oxygen atoms**.
2.  **Degree of Unsaturation (DoU):** The DoU represents the total number of rings and pi (π) bonds. The molecule has one C=O double bond and one C=C double bond. Therefore, the total DoU must be 1 + 1 = **2**.

**Step 3: Evaluate the Candidate Formulas**

Let's test the given options against our two criteria. The formula for DoU is `C + 1 - (H/2)`.

*   **A) C6H10O:**
    *   Oxygen Count: 1. **Fails.**
*   **B) C6H12O2:**
    *   Oxygen Count: 2. (Passes)
    *   DoU = 6 + 1 - (12/2) = 7 - 6 = 1. **Fails.** (This DoU only accounts for the C=O bond, not the C=C bond).
*   **C) C6H10O2:**
    *   Oxygen Count: 2. (Passes)
    *   DoU = 6 + 1 - (10/2) = 7 - 5 = 2. **Passes.** (This DoU perfectly accounts for both the C=O and C=C bonds).
*   **D) C6H12O:**
    *   Oxygen Count: 1. **Fails.**

The only formula that satisfies all experimental evidence is **C6H10O2**.

**Step 4: Final Answer Selection**

The correct molecular formula is C6H10O2. We must now match this to the correct letter from the question's list:
A) C6H10O
B) C6H12O2
C) C6H10O2
D) C6H12O

The formula C6H10O2 corresponds to option **C**.

**Careful Points:**
*   The broad peak at 3000 cm⁻¹ is key to identifying the O-H of a carboxylic acid, distinguishing it from a sharper alcohol peak.
*   The m/z = 45 fragment is a powerful confirmation of the carboxylic acid group.
*   Both the oxygen count and the Degree of Unsaturation must be used to differentiate between C6H12O2 and C6H10O2.
*   A common mistake is to perform the chemical analysis correctly but then select the wrong letter. It is crucial to carefully map the derived formula (C6H10O2) back to the specific options provided.

<<<C>>>
"""

# Run the check
result = check_answer(final_answer_text)
print(result)