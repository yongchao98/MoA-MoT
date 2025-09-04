import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by verifying the chosen chemical formula
    against constraints derived from the problem's spectroscopic data.
    """

    # --- 1. Define Constraints from the Question ---
    # Degree of Unsaturation (DoU): Aromatic ring (4) + C=C vinyl (1) + C=O ester (1) = 6
    required_dou = 6
    # Hydrogen Count from 1H NMR: 4 (aromatic) + 2 (vinyl) + 2 * 3 (methyl) = 12
    required_h_count = 12
    # Carbon Count from deduced fragments: C6 (ring) + C3 (propenyl) + C2 (methyl ester) = 11
    required_c_count = 11
    # Oxygen Count from ester group (-COO-)
    required_o_count = 2

    # --- 2. Define the LLM's Answer and the Options ---
    llm_answer_option = "B"
    options = {
        "A": "C12H12O2",
        "B": "C11H12O2",
        "C": "C12H14O2",
        "D": "C11H14O2"
    }
    
    chosen_formula = options.get(llm_answer_option)
    if not chosen_formula:
        return f"Invalid answer option '{llm_answer_option}'. It is not one of the provided choices."

    # --- 3. Helper Functions to Analyze a Formula ---
    def parse_formula(formula):
        """Parses a chemical formula string into a dictionary of atom counts."""
        counts = {'C': 0, 'H': 0, 'O': 0}
        parts = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        for element, count in parts:
            counts[element] = int(count) if count else 1
        return counts

    def calculate_dou(counts):
        """Calculates the Degree of Unsaturation."""
        # DoU = C + 1 - (H/2)
        return counts.get('C', 0) + 1 - (counts.get('H', 0) / 2)

    # --- 4. Verify the Chosen Formula ---
    atom_counts = parse_formula(chosen_formula)
    c = atom_counts.get('C', 0)
    h = atom_counts.get('H', 0)
    o = atom_counts.get('O', 0)
    dou = calculate_dou(atom_counts)

    # Check each constraint
    if o != required_o_count:
        return (f"Incorrect. The answer {chosen_formula} is wrong because it has {o} oxygen atoms, "
                f"but the presence of an ester group requires {required_o_count}.")

    if dou != required_dou:
        return (f"Incorrect. The answer {chosen_formula} is wrong because its Degree of Unsaturation (DoU) is {dou}, "
                f"but the structure described (aromatic ring + vinyl + ester) requires a DoU of {required_dou}.")

    if h != required_h_count:
        return (f"Incorrect. The answer {chosen_formula} is wrong because it has {h} hydrogen atoms, "
                f"but the 1H NMR data indicates a total of {required_h_count} hydrogens.")

    if c != required_c_count:
        return (f"Incorrect. The answer {chosen_formula} is wrong because it has {c} carbon atoms, "
                f"but the structure deduced from the NMR data requires {required_c_count} carbons.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_correctness()
print(result)