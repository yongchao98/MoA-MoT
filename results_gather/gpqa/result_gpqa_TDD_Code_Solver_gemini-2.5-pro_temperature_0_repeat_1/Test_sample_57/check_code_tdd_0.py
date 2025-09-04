import re

def check_answer():
    """
    Checks the correctness of the given answer to the physics question.
    """
    question = "Which of the following physical theories never requires regularization at high energies?"
    options = {
        "A": "Superstring Theory",
        "B": "Quantum Chromodynamics",
        "C": "Classical Electrodynamics",
        "D": "Quantum Electrodynamics"
    }

    # Physics knowledge base: True if the theory requires high-energy regularization.
    regularization_needed = {
        "Superstring Theory": False,
        "Quantum Chromodynamics": True,
        "Classical Electrodynamics": True, # Due to the self-energy problem of a point charge.
        "Quantum Electrodynamics": True
    }

    # The provided response from the LLM
    llm_response = "<<<A>>>"

    # Extract the answer from the response
    match = re.search(r'<<<([A-Z])>>>', llm_response)
    if not match:
        return "Incorrect: The answer format is invalid. It should be like '<<<A>>>'."

    llm_answer_key = match.group(1)

    # Check if the answer key is a valid option
    if llm_answer_key not in options:
        return f"Incorrect: The answer '{llm_answer_key}' is not a valid option. Options are A, B, C, D."

    # Get the theory corresponding to the LLM's answer
    chosen_theory = options[llm_answer_key]

    # Find the correct answer key based on our knowledge base
    correct_key = None
    for key, theory in options.items():
        if not regularization_needed[theory]:
            correct_key = key
            break
    
    correct_theory = options[correct_key]

    # Check if the LLM's answer is correct
    if llm_answer_key == correct_key:
        return "Correct"
    else:
        return (f"Incorrect: The chosen answer '{llm_answer_key}' corresponds to '{chosen_theory}'. "
                f"This theory requires high-energy regularization. "
                f"The correct answer is '{correct_key}' ({correct_theory}), which is conjectured to be UV-finite and does not require regularization.")

# Run the check
result = check_answer()
print(result)