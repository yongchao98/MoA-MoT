import re

def check_chemistry_answer():
    """
    This function verifies the correct answer for the given chemistry problem
    by applying constraints derived from the spectroscopic data.
    """
    # --- Step 1: Define the options from the question ---
    options = {
        'A': 'C6H10O2',
        'B': 'C6H12O',
        'C': 'C6H10O',
        'D': 'C6H12O2'
    }

    # --- Step 2: Define constraints based on the problem description ---
    # Constraint 1: Carboxylic acid (-COOH) requires 2 Oxygen atoms.
    required_oxygens = 2
    # Constraint 2: Carboxylic acid (C=O) + Alkene (C=C) requires 2 degrees of unsaturation.
    required_dou = 2

    # --- Step 3: Helper functions to parse formula and calculate DoU ---
    def parse_formula(formula_str):
        """Parses a chemical formula string into a dictionary of atom counts."""
        atoms = {'C': 0, 'H': 0, 'O': 0}
        # Find all occurrences of an element followed by an optional number
        for element, count in re.findall(r'([A-Z])(\d*)', formula_str):
            atoms[element] = int(count) if count else 1
        return atoms

    def calculate_dou(c_atoms, h_atoms):
        """Calculates the Degree of Unsaturation for a C, H, O compound."""
        # Formula: DoU = C + 1 - (H/2)
        return c_atoms + 1 - (h_atoms / 2)

    # --- Step 4: Iterate through options and check against constraints ---
    correct_options = []
    print("--- Verifying Options Against Spectroscopic Data ---")
    for option_letter, formula in options.items():
        atoms = parse_formula(formula)
        num_oxygens = atoms.get('O', 0)
        dou = calculate_dou(atoms.get('C', 0), atoms.get('H', 0))

        # Check if both constraints are met
        oxygen_match = (num_oxygens == required_oxygens)
        dou_match = (dou == required_dou)

        print(f"\nChecking Option {option_letter}: {formula}")
        print(f"  - Oxygen count: {num_oxygens} (Required: {required_oxygens}) -> {'Pass' if oxygen_match else 'Fail'}")
        print(f"  - Degree of Unsaturation: {dou} (Required: {required_dou}) -> {'Pass' if dou_match else 'Fail'}")

        if oxygen_match and dou_match:
            correct_options.append(option_letter)
            print("  Result: This option satisfies all constraints.")
        else:
            print("  Result: This option is incorrect.")


    # --- Step 5: Report the final correct answer ---
    if len(correct_options) == 1:
        final_answer = correct_options[0]
        print(f"\nConclusion: The only correct answer is Option {final_answer} ({options[final_answer]}).")
        return "Correct" # Simulating the check for the correct answer
    elif len(correct_options) == 0:
        print("\nConclusion: No option satisfies all the constraints.")
        return "Incorrect: No valid option found."
    else:
        print(f"\nConclusion: Multiple options satisfy the constraints: {correct_options}")
        return "Incorrect: Ambiguous result, multiple options are valid."

# To run the check, you would call the function.
# The output of this function confirms the analysis.
check_chemistry_answer()