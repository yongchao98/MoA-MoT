def check_genetics_question_answer(llm_answer: str) -> str:
    """
    Checks the correctness of the answer to a molecular genetics question.

    The function models the biological constraints as a set of rules and
    evaluates each option against them.

    Args:
        llm_answer: The final answer provided by the LLM (e.g., 'A', 'B', 'C', 'D').

    Returns:
        A string indicating "Correct" or the reason for being incorrect.
    """

    # Define the options with their key characteristics based on the question's context.
    options = {
        'A': {
            'description': 'loss of protein dimerization and wild-type phenotype',
            'phenotype': 'wild-type',
            'interaction': 'lost'
        },
        'B': {
            'description': 'protein aggregation and loss-of-function phenotype',
            'phenotype': 'loss-of-function',
            'interaction': 'present' # Aggregation implies interaction
        },
        'C': {
            'description': 'change of protein conformation and gain-of-function phenotype',
            'phenotype': 'gain-of-function',
            'interaction': 'present'
        },
        'D': {
            'description': 'protein degradation and loss-of-function of the wild-type allele',
            'phenotype': 'loss-of-function of the wild-type allele', # A specific and defining type of loss-of-function
            'interaction': 'present' # Degradation of the heterodimer implies interaction
        }
    }

    # --- Rule-based evaluation ---
    reasons_for_rejection = {}
    valid_options = []

    for label, props in options.items():
        is_valid = True
        rejection_log = []

        # Rule 1: Check the phenotype. Must be loss-of-function.
        if 'gain-of-function' in props['phenotype'] or 'wild-type' in props['phenotype']:
            is_valid = False
            rejection_log.append(f"A dominant-negative mutation causes a loss-of-function, not a '{props['phenotype']}'.")

        # Rule 2: Check the interaction mechanism. Must allow for interference.
        if props['interaction'] == 'lost':
            is_valid = False
            rejection_log.append("A dominant-negative effect requires the mutant protein to interact with the wild-type protein. A complete loss of dimerization would prevent this interference.")

        if is_valid:
            valid_options.append(label)
        else:
            reasons_for_rejection[label] = " ".join(rejection_log)

    # --- Nuanced comparison of remaining valid options (B and D) ---
    # Rule 3: Check for the most precise and fundamental description.
    best_option = None
    if 'B' in valid_options and 'D' in valid_options:
        # Option D's "loss-of-function of the wild-type allele" is the textbook definition
        # of a dominant-negative effect. This is more precise than the general
        # "loss-of-function phenotype" in option B.
        best_option = 'D'
    elif len(valid_options) == 1:
        best_option = valid_options[0]

    # --- Final Verdict ---
    if llm_answer == best_option:
        return "Correct"
    else:
        if llm_answer in reasons_for_rejection:
            return f"Incorrect. The provided answer '{llm_answer}' violates fundamental principles: {reasons_for_rejection[llm_answer]}"
        elif best_option:
            return f"Incorrect. The provided answer is '{llm_answer}', but the best answer is '{best_option}'. Option D is superior because 'loss-of-function of the wild-type allele' is the precise definition of a dominant-negative effect, making it a more accurate and fundamental description than the other plausible options."
        else:
            return f"Incorrect. The provided answer '{llm_answer}' is wrong, and the logic could not determine a single best answer from the options."

# The final answer from the analysis is 'D'.
provided_answer = 'D'
result = check_genetics_question_answer(provided_answer)
print(result)