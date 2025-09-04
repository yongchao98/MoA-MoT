import re

def check_correctness_of_answer(llm_response):
    """
    Checks if the provided LLM answer to the physics question is correct.

    The question is:
    Which of the following physical theories never requires regularization at high energies?
    A) Superstring Theory
    B) Quantum Chromodynamics
    C) Classical Electrodynamics
    D) Quantum Electrodynamics

    The correct answer is A) Superstring Theory.
    - Quantum Field Theories (QFTs) like Quantum Chromodynamics (B) and Quantum Electrodynamics (D)
      suffer from ultraviolet (high-energy) divergences in calculations and require regularization
      and renormalization to yield finite, physical predictions.
    - Classical Electrodynamics (C) has its own divergence problem related to the infinite self-energy
      of a point charge, which also requires a form of regularization (e.g., assuming a finite particle size).
    - Superstring Theory (A) is unique in this list as it is conjectured to be "UV-finite". By replacing
      point-like particles with extended one-dimensional strings, it is believed to resolve the divergences
      that plague QFTs, thus not requiring regularization.
    """
    
    # Define the mapping of options to theories and the correct option.
    options = {
        "A": "Superstring Theory",
        "B": "Quantum Chromodynamics",
        "C": "Classical Electrodynamics",
        "D": "Quantum Electrodynamics"
    }
    correct_option = "A"
    correct_theory = options[correct_option]

    # Extract the final answer from the LLM's response using regex.
    # The expected format is <<<X>>> where X is the option letter.
    match = re.search(r'<<<([A-D])>>>', llm_response)
    
    if not match:
        return "The answer format is incorrect. The final answer should be enclosed in '<<<' and '>>>', for example, '<<<A>>>'."

    llm_selected_option = match.group(1)

    # Check if the selected option is the correct one.
    if llm_selected_option == correct_option:
        return "Correct"
    else:
        llm_selected_theory = options.get(llm_selected_option, "an invalid option")
        reason = (
            f"The answer is incorrect. The LLM chose option {llm_selected_option} ({llm_selected_theory}).\n"
            f"The correct answer is option {correct_option} ({correct_theory}).\n"
            "Reasoning: Superstring Theory is the only theory listed that is conjectured to be UV-finite, "
            "meaning it does not require regularization for high-energy calculations. "
            "Quantum Chromodynamics, Quantum Electrodynamics, and Classical Electrodynamics all contain "
            "divergences that must be handled with regularization techniques."
        )
        return reason

# The user-provided response from the other LLM
llm_answer_from_user = """
The checker reported a syntax error. This was likely caused by special characters in the explanatory text of my previous response. I will now provide the code solution again, ensuring the surrounding text is simple to avoid parsing issues. The underlying code logic, which correctly identifies Superstring Theory, remains unchanged.
