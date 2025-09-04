import re

def check_answer(llm_response):
    """
    Checks the correctness of the LLM's answer for the SMEFT symmetry question.

    The function verifies that the chosen option corresponds to the set of symmetries
    that must be respected by all operators in the Standard Model Effective Field Theory (SMEFT).
    """
    # --- Define the ground truth based on physics principles ---
    # 1. Lorentz Symmetry: Required. SMEFT is a relativistic QFT.
    # 2. Poincaré Symmetry: Required. This is the full spacetime symmetry group for a relativistic QFT.
    # 3. CP Symmetry: Not required. The SM itself violates CP, and SMEFT is used to parameterize new sources of CP violation.
    # 4. CPT Symmetry: Required. A consequence of the CPT theorem for any local, Lorentz-invariant QFT.
    correct_symmetries = {1, 2, 4}
    
    # --- Define the options as presented in the question ---
    options = {
        'A': {1, 3, 4},
        'B': {1, 2},
        'C': {3, 4},
        'D': {1, 2, 4}
    }
    
    # --- Extract the final answer from the LLM's response ---
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Invalid answer format. The final answer should be in the format <<<X>>> where X is A, B, C, or D."
        
    answer_key = match.group(1)
    
    # --- Compare the LLM's choice with the ground truth ---
    selected_symmetries = options.get(answer_key)
    
    if selected_symmetries == correct_symmetries:
        return "Correct"
    else:
        # --- Generate a detailed reason for the incorrect answer ---
        reason = f"The answer '{answer_key}' is incorrect. "
        
        missing_symmetries = correct_symmetries - selected_symmetries
        if missing_symmetries:
            reason += f"It fails to include the required symmetry/symmetries: {sorted(list(missing_symmetries))}. "
            
        extra_symmetries = selected_symmetries - correct_symmetries
        if extra_symmetries:
            reason += f"It incorrectly includes symmetry {extra_symmetries.pop()}, which is not a required symmetry for all SMEFT operators. "
            
        reason += f"The correct set of required symmetries is {sorted(list(correct_symmetries))}, corresponding to Lorentz, Poincaré, and CPT symmetry."
        return reason

# The final answer provided by the LLM
llm_final_answer = "<<<D>>>"

# Check the correctness of the answer
result = check_answer(llm_final_answer)
print(result)