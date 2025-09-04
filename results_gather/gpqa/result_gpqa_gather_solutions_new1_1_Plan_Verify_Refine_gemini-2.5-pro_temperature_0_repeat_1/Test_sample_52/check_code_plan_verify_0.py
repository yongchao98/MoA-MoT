import re

def check_answer():
    """
    Checks the correctness of the LLM's answer by deducing the chemical formula
    from the spectroscopic data provided in the question.
    """

    # 1. Define constraints from the spectroscopic data.
    # - Aromatic ring (1 ring, 3 double bonds) -> DoU = 4
    # - Ester group (-COO-) has a C=O double bond -> DoU = 1
    # - Vinyl-H signals indicate a C=C double bond -> DoU = 1
    required_dou = 4 + 1 + 1

    # The "no -CH2- groups" is a crucial structural constraint.
    # The fragments are: C6H4 (aromatic), C3H5 (propenyl), and C2H3O2 (methyl ester).
    # This leads to a base formula of C11H12O2.

    # 2. Define the options from the question.
    options = {
        "A": "C11H12O2",
        "B": "C12H14O2",
        "C": "C12H12O2",
        "D": "C11H14O2",
    }

    # 3. The LLM's answer to be checked.
    llm_answer_letter = "A"
    llm_formula = options.get(llm_answer_letter)

    if not llm_formula:
        return f"Invalid answer format. The answer '{llm_answer_letter}' does not correspond to any option."

    # 4. Helper function to parse formula and calculate DoU.
    def get_formula_info(formula_str):
        atoms = {'C': 0, 'H': 0, 'O': 0, 'N': 0}
        for element, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula_str):
            atoms[element] = int(count) if count else 1
        
        # Formula for Degree of Unsaturation (DoU)
        dou = atoms['C'] - (atoms['H'] / 2) + (atoms['N'] / 2) + 1
        return atoms, dou

    # 5. Systematically check all options to find the one that fits all constraints.
    correct_candidates = []
    reasons_for_failure = {}

    for letter, formula in options.items():
        atoms, dou = get_formula_info(formula)
        
        # Constraint 1: Check Degree of Unsaturation
        if dou != required_dou:
            reasons_for_failure[letter] = f"Formula {formula} has a Degree of Unsaturation (DoU) of {dou}, but the structure requires a DoU of {required_dou}."
            continue

        # Constraint 2: Check for implied -CH2- groups.
        # The base structure derived from fragments is C11H12O2.
        # We compare other formulas with the correct DoU to this base.
        if formula == "C12H14O2":
            # This formula has an extra CH2 compared to the base C11H12O2.
            reasons_for_failure[letter] = f"Formula {formula} has the correct DoU, but it implies an extra -CH2- group compared to the structure derived from the NMR fragments (C11H12O2), which violates the 'no -CH2- groups' constraint."
            continue
            
        # If a formula passes all checks, it's a correct candidate.
        correct_candidates.append(letter)

    # 6. Final validation and result.
    if len(correct_candidates) != 1:
        return f"Analysis failed to isolate a single correct answer. Candidates found: {correct_candidates}. Reasons for failure: {reasons_for_failure}"

    correct_letter = correct_candidates[0]

    if llm_answer_letter == correct_letter:
        return "Correct"
    else:
        failure_reason = reasons_for_failure.get(llm_answer_letter, 
            f"The answer {llm_formula} does not satisfy all constraints.")
        return f"Incorrect. The correct answer is {correct_letter} ({options[correct_letter]}). The provided answer {llm_answer_letter} ({llm_formula}) is wrong. Reason: {failure_reason}"

# Execute the check
result = check_answer()
print(result)