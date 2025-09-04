def check_polymerization_answer():
    """
    Checks the correctness of the final answer for the ethylene polymerization question.
    The "ground truth" is synthesized from the consensus and analysis of the provided LLM answers.
    """
    # The final answer provided by the user's system to be checked.
    final_answer = "<<<C>>>"
    
    # --- Ground Truth Database ---
    # This dictionary represents the established facts based on a consensus of chemical literature
    # and the provided rationales from various LLMs.
    statement_facts = {
        'A': {
            "is_correct": True,
            "is_best_answer": False,
            "reason": "Statement A ('Certain noble metal catalysts can be used but are too expensive.') is a plausible statement. However, these catalysts (e.g., Palladium-based) typically produce hyperbranched structures via a 'chain-walking' mechanism, which does not fit the question's constraint of 'regular branches' that are formed by incorporating a specific comonomer. Therefore, it is not the best answer."
        },
        'B': {
            "is_correct": False,
            "is_best_answer": False,
            "reason": "Statement B ('Aluminum-based activators do not work for the essential additional reaction step.') is factually incorrect. The oligomerization step requires an activator, and aluminum-based compounds like methylaluminoxane (MAO) are the most common and essential activators for the relevant chromium-based catalysts."
        },
        'C': {
            "is_correct": True,
            "is_best_answer": True,
            "reason": "Statement C ('One can use a catalyst of a group VIa transition metal in combination with specific activators.') is the most accurate and fundamental chemical statement. Chromium (Cr), a Group VIa metal, is the key component in state-of-the-art industrial catalysts that selectively convert ethylene to 1-hexene or 1-octene, which is the essential step for creating regular branches."
        },
        'D': {
            "is_correct": False,
            "is_best_answer": False,
            "reason": "Statement D ('Such combined systems are already implemented on an industrial scale in the US.') is technically incorrect under a strict interpretation. While the component technologies (e.g., 1-hexene production) are industrialized in the US, a true single-reactor 'combined system' is not in widespread industrial use due to differing optimal conditions for the two catalysts. The statement is an oversimplification."
        }
    }

    # --- Verification Logic ---
    try:
        chosen_option = final_answer.strip().replace('<', '').replace('>', '')
        if chosen_option not in statement_facts:
            return f"Error: Invalid option '{chosen_option}'. The option must be one of A, B, C, D."
    except Exception:
        return "Error: Invalid answer format. Expected format like '<<<C>>>'."

    # Check if the chosen option is considered the single best answer.
    if statement_facts[chosen_option]['is_best_answer']:
        return "Correct"
    else:
        # If the chosen answer is not the best, explain why.
        # First, state why the chosen option is wrong.
        reason_for_rejection = f"The provided answer '{chosen_option}' is incorrect. Reason: {statement_facts[chosen_option]['reason']}"
        
        # Then, identify the correct answer and explain why it's correct.
        correct_option = None
        for option, details in statement_facts.items():
            if details['is_best_answer']:
                correct_option = option
                break
        
        if correct_option:
            reason_for_correctness = f"The correct answer should be '{correct_option}'. Reason: {statement_facts[correct_option]['reason']}"
            return f"{reason_for_rejection}\n\n{reason_for_correctness}"
        else:
            # This case should not be reached if the facts are set up correctly with one best answer.
            return f"{reason_for_rejection}\n\nError: No single best answer could be identified from the facts."

# Execute the check and print the result.
result = check_polymerization_answer()
print(result)