import re

def check_smeft_symmetry_answer():
    """
    Checks the correctness of the answer about SMEFT symmetries.

    The function encodes the physical principles of the Standard Model Effective Field Theory (SMEFT)
    to verify the provided answer.
    """
    # Define the question's options and the final answer from the LLM
    question_options = {
        'A': {1, 2, 4},
        'B': {3, 4},
        'C': {1, 2},
        'D': {1, 3, 4}
    }
    
    symmetries_map = {
        1: "Lorentz Symmetry",
        2: "Poincare symmetry",
        3: "CP symmetry",
        4: "CPT symmetry"
    }

    final_answer_text = "<<<A>>>"

    # --- Ground Truth based on Physics Principles ---
    # 1. Lorentz Symmetry: Required. SMEFT is a relativistic QFT.
    # 2. Poincare Symmetry: Required. This is the full spacetime symmetry group of special relativity.
    # 3. CP Symmetry: Not required. The SM itself violates CP, and SMEFT is designed to parameterize new sources of CP violation.
    # 4. CPT Symmetry: Required. This is a consequence of the CPT theorem for local, Lorentz-invariant QFTs.
    correct_symmetries = {1, 2, 4}

    # Extract the letter from the final answer
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Invalid answer format. The answer should be in the format <<<X>>> where X is A, B, C, or D."
    
    chosen_letter = match.group(1)
    
    if chosen_letter not in question_options:
        return f"The chosen answer letter '{chosen_letter}' is not a valid option."

    chosen_symmetries = question_options[chosen_letter]

    # --- Verification Logic ---
    if chosen_symmetries == correct_symmetries:
        return "Correct"
    else:
        # Analyze the error
        missing_symmetries = correct_symmetries - chosen_symmetries
        extra_symmetries = chosen_symmetries - correct_symmetries
        
        error_messages = []
        if missing_symmetries:
            missing_names = [f"{s} ({symmetries_map[s]})" for s in sorted(list(missing_symmetries))]
            error_messages.append(f"The answer is incorrect because it is missing the required symmetry/symmetries: {', '.join(missing_names)}.")
        
        if extra_symmetries:
            extra_names = [f"{s} ({symmetries_map[s]})" for s in sorted(list(extra_symmetries))]
            error_messages.append(f"The answer is incorrect because it includes the symmetry/symmetries: {', '.join(extra_names)}, which are not required to be respected by all SMEFT operators.")
            
        return "\n".join(error_messages)

# Execute the check and print the result
result = check_smeft_symmetry_answer()
print(result)