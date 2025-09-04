import re

def check_answer_correctness(llm_response: str) -> str:
    """
    Checks the correctness of the answer to the physics question.

    Question: Which of the following physical theories never requires regularization at high energies?
        A) Superstring Theory
        B) Quantum Chromodynamics
        C) Classical Electrodynamics
        D) Quantum Electrodynamics

    Args:
        llm_response: The string containing the LLM's full response, including the final answer.

    Returns:
        A string indicating "Correct" or the reason for the error.
    """
    # Based on physics principles, Superstring Theory is the correct answer.
    # It is conjectured to be UV-finite, unlike the other options which all suffer
    # from high-energy / short-distance divergences.
    correct_answer_key = "A"
    
    # Use regex to find the answer in the format <<<X>>>
    match = re.search(r'<<<([A-D])>>>', llm_response)

    if not match:
        return "The answer format is incorrect. The final answer must be in the format '<<<X>>>' where X is A, B, C, or D."

    provided_answer_key = match.group(1)

    if provided_answer_key == correct_answer_key:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer was '{provided_answer_key}', but the correct answer is '{correct_answer_key}'.\n"
            "Reasoning: Superstring Theory (A) is conjectured to be UV-finite and does not require regularization. "
            "Quantum Chromodynamics (B) and Quantum Electrodynamics (D) are quantum field theories that famously require regularization. "
            "Classical Electrodynamics (C) also has a divergence problem in the self-energy of a point charge, which is a high-energy issue."
        )
        return reason

# The user-provided response to be checked
llm_response_to_check = """
The checker reported a `SyntaxError: unterminated triple-quoted string literal`. This error occurs when the code I provide contains a `"""` sequence, which prematurely closes the string literal in the checker's validation script. My previous solution used a multi-line docstring enclosed in `"""`, which caused this issue.

To fix this, I will rewrite the code, replacing the multi-line docstring with single-line `#` comments. The function logic, tests, and final answer remain the same, as they were correct. This change only addresses the formatting to avoid the syntax error.
