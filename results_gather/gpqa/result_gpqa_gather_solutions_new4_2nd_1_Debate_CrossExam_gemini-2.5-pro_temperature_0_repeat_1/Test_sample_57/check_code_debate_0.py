import re

def check_physics_theory_answer(question_options, llm_answer_text):
    """
    Checks the correctness of the answer to the physics theory question.

    Args:
        question_options (dict): A dictionary mapping option letters to theory names.
        llm_answer_text (str): The text of the LLM's final answer.

    Returns:
        str: "Correct" if the answer is correct, otherwise a reason for the error.
    """

    # Knowledge base about which theories require regularization.
    # True = Requires regularization
    # False = Does not require regularization
    theory_properties = {
        "Quantum Electrodynamics": {
            "requires_regularization": True,
            "reason": "it is a quantum field theory with UV divergences in loop diagrams."
        },
        "Classical Electrodynamics": {
            "requires_regularization": True,
            "reason": "it has the infinite self-energy problem for point charges, a short-distance divergence."
        },
        "Superstring Theory": {
            "requires_regularization": False,
            "reason": "it is believed to be UV-finite because the extended nature of strings smooths out interactions at short distances."
        },
        "Quantum Chromodynamics": {
            "requires_regularization": True,
            "reason": "it is a quantum field theory with UV divergences in loop diagrams."
        }
    }

    # Find the correct answer from the knowledge base
    correct_theory_name = None
    for theory, properties in theory_properties.items():
        if not properties["requires_regularization"]:
            correct_theory_name = theory
            break
    
    if not correct_theory_name:
        return "Error in checker: No correct theory found in the knowledge base."

    # Find the correct option letter from the question's mapping
    correct_option_letter = None
    for letter, theory in question_options.items():
        if theory == correct_theory_name:
            correct_option_letter = letter
            break

    if not correct_option_letter:
        return f"Error in checker: The correct theory '{correct_theory_name}' is not in the provided options."

    # --- Parse the LLM's answer ---
    # The final answer is usually on the last line, e.g., "C) Superstring Theory"
    # We also handle the <<<...>>> format if present.
    final_answer_line = ""
    if "<<<" in llm_answer_text:
        match = re.search(r'<<<(.+?)>>>', llm_answer_text, re.DOTALL)
        if match:
            final_answer_line = match.group(1).strip()
        else: # Fallback to last line if <<<>>> is malformed
            final_answer_line = llm_answer_text.strip().split('\n')[-1].strip()
    else:
        final_answer_line = llm_answer_text.strip().split('\n')[-1].strip()

    # Extract the letter and theory name from the parsed line
    # e.g., "C) Superstring Theory" -> ("C", "Superstring Theory")
    match = re.match(r'([A-D])\)\s*(.*)', final_answer_line)
    if not match:
        return f"Incorrect. The answer format is not recognized. Could not parse '{final_answer_line}'."

    provided_option_letter = match.group(1).strip()
    provided_theory_name = match.group(2).strip()

    # --- Perform Checks ---
    # 1. Check if the provided theory name is correct
    if provided_theory_name != correct_theory_name:
        reason_for_incorrectness = theory_properties.get(provided_theory_name, {}).get("reason", "it is not the correct choice.")
        return (f"Incorrect. The answer identifies '{provided_theory_name}' as the correct theory, but it is wrong because {reason_for_incorrectness} "
                f"The correct theory is '{correct_theory_name}'.")

    # 2. Check if the provided option letter is correct
    if provided_option_letter != correct_option_letter:
        return (f"Incorrect. While '{provided_theory_name}' is the correct theory, the answer incorrectly labels it as option '{provided_option_letter}'. "
                f"According to the question, it should be option '{correct_option_letter}'.")

    # If both checks pass, the answer is correct.
    return "Correct"

# --- Main execution ---
# Define the options as they were presented in the original question
# This is crucial because different LLMs might have reordered them in their analysis.
original_question_options = {
    "A": "Quantum Electrodynamics",
    "B": "Classical Electrodynamics",
    "C": "Superstring Theory",
    "D": "Quantum Chromodynamics"
}

# The final answer provided by the LLM to be checked.
llm_final_answer = """
Based on a step-by-step analysis of the provided answers and established physical principles, the correct answer is determined as follows:

**1. Understand the Core Concept: Regularization at High Energies**
The question asks which physical theory is free from "ultraviolet (UV) divergences," which are infinities that appear in calculations involving very high energies or, equivalently, very short distances. "Regularization" is the mathematical technique used to manage these infinities. A theory that "never requires regularization" must be inherently free from these types of infinities.

**2. Analyze Each Option based on Consensus Reasoning**
An analysis of the reasoning within the provided candidate answers reveals a strong, near-unanimous consensus on the physics of each theory, even though their final lettered choices are contradictory. This consensus is confirmed by external search results.

*   **A) Quantum Electrodynamics (QED):** This is the quantum theory of electromagnetism. It is the classic example of a theory that *requires* regularization. Its calculations produce infinities that must be managed through regularization and renormalization to yield famously precise predictions.
*   **B) Classical Electrodynamics:** While not a quantum theory, it has its own high-energy problem: the self-energy of a point-like charge is infinite. This divergence shows the theory breaks down at short distances and requires a form of regularization to be resolved.
*   **D) Quantum Chromodynamics (QCD):** This is the quantum theory of the strong force. As a quantum field theory similar to QED, it also suffers from UV divergences in its calculations. Therefore, QCD absolutely requires regularization to be a predictive theory.
*   **C) Superstring Theory:** A primary motivation for superstring theory was to resolve the problem of infinities. By modeling fundamental particles as one-dimensional "strings" rather than zero-dimensional points, the theory "smears out" interactions. This intrinsic feature naturally eliminates the short-distance/high-energy infinities found in the other theories. For this reason, superstring theory is believed to be "UV-finite" and does not require regularization for high-energy interactions.

**3. Conclusion**
The reasoning presented in the vast majority of the candidate answers, along with established physics principles, consistently identifies that QED, QCD, and Classical Electrodynamics all require regularization. Superstring Theory is the only option specifically formulated to be free of these high-energy divergences. Therefore, it is the correct answer.

C) Superstring Theory
"""

# Run the check and print the result
result = check_physics_theory_answer(original_question_options, llm_final_answer)
print(result)