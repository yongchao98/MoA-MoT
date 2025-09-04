import re

def check_physics_theory_answer():
    """
    Checks the correctness of the answer to the physics question.

    The question asks which theory never requires regularization at high energies.
    This check is based on established principles of theoretical physics.
    """

    # The options as defined in the final provided answer being checked
    options_map = {
        "A": "Superstring Theory",
        "B": "Classical Electrodynamics",
        "C": "Quantum Electrodynamics",
        "D": "Quantum Chromodynamics"
    }

    # The final answer provided by the LLM to be checked
    llm_answer = "<<<A>>>"

    # Knowledge base: A dictionary mapping each theory to whether it requires regularization
    # and the reason why. True means it requires regularization.
    theory_properties = {
        "Superstring Theory": {
            "requires_regularization": False,
            "reason": "It is believed to be UV-finite. The extended, non-point-like nature of strings 'smears out' interactions, naturally avoiding the high-energy (ultraviolet) divergences found in point-particle theories."
        },
        "Classical Electrodynamics": {
            "requires_regularization": True,
            "reason": "It suffers from the infinite self-energy problem for a point-like charge, a high-energy/short-distance divergence that requires a form of regularization (e.g., a cutoff) to resolve."
        },
        "Quantum Electrodynamics": {
            "requires_regularization": True,
            "reason": "It is a quantum field theory that famously requires regularization and renormalization to handle the ultraviolet (UV) divergences that appear in calculations involving loop diagrams."
        },
        "Quantum Chromodynamics": {
            "requires_regularization": True,
            "reason": "Like QED, it is a quantum field theory whose calculations produce UV divergences, making regularization and renormalization essential procedures to obtain finite, physical predictions."
        }
    }

    # Determine the correct option based on the knowledge base
    correct_option_letter = None
    for letter, theory_name in options_map.items():
        if not theory_properties[theory_name]["requires_regularization"]:
            correct_option_letter = letter
            break
    
    if correct_option_letter is None:
        # This case should not be reached with the current knowledge base
        return "Error in checking logic: No correct answer found in the knowledge base."

    # Extract the letter from the LLM's answer string
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return f"Invalid answer format. Expected '<<<X>>>' where X is A, B, C, or D. Got: {llm_answer}"
    
    llm_choice_letter = match.group(1)

    # Compare the LLM's choice with the correct answer
    if llm_choice_letter == correct_option_letter:
        return "Correct"
    else:
        chosen_theory = options_map.get(llm_choice_letter, "Invalid Option")
        correct_theory = options_map.get(correct_option_letter, "N/A")

        error_message = (
            f"Incorrect. The provided answer is {llm_choice_letter} ({chosen_theory}), but the correct answer is {correct_option_letter} ({correct_theory}).\n\n"
            f"Reasoning:\n"
            f"The question asks for the theory that does NOT require regularization at high energies.\n"
            f"- The chosen answer, {chosen_theory}, is incorrect because: {theory_properties[chosen_theory]['reason']}\n"
            f"- The correct answer, {correct_theory}, is correct because: {theory_properties[correct_theory]['reason']}"
        )
        return error_message

# Run the check and print the result
result = check_physics_theory_answer()
print(result)