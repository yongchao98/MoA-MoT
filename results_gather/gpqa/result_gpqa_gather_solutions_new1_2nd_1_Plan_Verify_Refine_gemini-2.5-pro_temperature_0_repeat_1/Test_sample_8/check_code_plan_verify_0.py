import re

def check_molecular_symmetry_answer():
    """
    Checks the correctness of the LLM's answer to the molecular symmetry question
    by comparing it against established chemical knowledge.
    """
    # Ground truth data: A dictionary mapping molecule names to their point groups.
    # For flexible molecules, the point group of the most stable, lowest-energy conformation is used.
    molecule_point_groups = {
        "quinuclidine": "C3v",
        "triphenyleno[1,2-c:5,6-c':9,10-c'']trifuran-1,3,6,8,11,13-hexaone": "C3h",
        "triisopropyl borate": "C3",  # The most stable conformation is a C3 propeller to minimize steric hindrance, not C3h.
        "benzo[1,2-c:3,4-c':5,6-c'']trifuran-1,3,4,6,7,9-hexaone": "D3h"
    }

    # The options as defined in the original question prompt.
    question_options = {
        "A": "quinuclidine",
        "B": "triphenyleno[1,2-c:5,6-c':9,10-c'']trifuran-1,3,6,8,11,13-hexaone",
        "C": "triisopropyl borate",
        "D": "benzo[1,2-c:3,4-c':5,6-c'']trifuran-1,3,4,6,7,9-hexaone"
    }

    # The final answer provided by the LLM to be checked.
    llm_answer_text = "<<<B>>>"

    # --- Verification Logic ---

    # 1. Extract the letter from the answer format.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Invalid answer format. Could not parse a letter from '{llm_answer_text}'."
    
    provided_answer_letter = match.group(1)

    # 2. Identify the molecule selected by the answer.
    selected_molecule_name = question_options.get(provided_answer_letter)
    if not selected_molecule_name:
        return f"Invalid answer. The letter '{provided_answer_letter}' is not one of the question options (A, B, C, D)."

    # 3. Find the correct molecule with C3h symmetry from our ground truth.
    correct_molecule = None
    correct_answer_letter = None
    for letter, molecule in question_options.items():
        if molecule_point_groups.get(molecule) == "C3h":
            correct_molecule = molecule
            correct_answer_letter = letter
            break
    
    if correct_molecule is None:
        return "Checker Error: No molecule with C3h symmetry found in the ground truth data for the given options."

    # 4. Compare the selected molecule with the correct molecule.
    if selected_molecule_name == correct_molecule:
        # The answer is correct. As a final sanity check, ensure no other options were also correct.
        other_correct_options = []
        for letter, molecule in question_options.items():
            if molecule != correct_molecule and molecule_point_groups.get(molecule) == "C3h":
                other_correct_options.append(letter)
        
        if not other_correct_options:
            return "Correct"
        else:
            # This case means the question had multiple correct answers.
            return f"The selected answer {provided_answer_letter} is correct, but the question is ambiguous as option(s) {', '.join(other_correct_options)} are also correct."
    else:
        # The answer is incorrect. Provide the reason.
        actual_point_group = molecule_point_groups.get(selected_molecule_name, "Unknown")
        reason = (f"Incorrect. The selected molecule (Option {provided_answer_letter}: {selected_molecule_name}) "
                  f"does not have C3h symmetry. Its point group is {actual_point_group}. "
                  f"The correct answer is Option {correct_answer_letter} ({correct_molecule}), which has C3h symmetry.")
        return reason

# Execute the check and print the result.
result = check_molecular_symmetry_answer()
print(result)