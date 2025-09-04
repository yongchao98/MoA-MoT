def check_biology_question_answer(llm_answer_choice: str):
    """
    Checks the correctness of an answer to a conceptual biology question
    by encoding the biological principles as a set of logical rules.

    Args:
        llm_answer_choice: The letter choice of the provided answer (e.g., "A", "B", "C", "D").

    Returns:
        A string indicating "Correct" or the reason for the incorrectness.
    """

    # --- Logical Model ---
    # 1. Define the scenario based on the question.
    # The protein must form a dimer to function.
    # A heterozygote produces both Wild-Type (WT) and Mutant (Y) protein subunits.
    # The mutation (Y) is in the dimerization domain.
    # The mutation is "dominant-negative".

    # 2. Deduce the molecular mechanism from the definitions.
    # "Dominant-negative" implies the mutant Y protein must interfere with the WT protein.
    # Since the mutation is in the dimerization domain and the protein functions as a dimer,
    # the interference must occur during dimerization.
    # This means Y can still bind to WT, forming a non-functional WT-Y heterodimer.
    # This process "poisons" or sequesters the WT protein, preventing it from forming
    # functional WT-WT homodimers.
    # The result is a significant reduction in overall protein function, i.e., a loss-of-function phenotype.

    # 3. Evaluate each option against the deduced mechanism.
    analysis = {
        "A": {
            "is_correct": False,
            "reason": "This is incorrect. For a dominant-negative effect to occur, the mutant protein must be stable enough to interact with and inhibit the wild-type protein. If the mutant protein were simply degraded, it could not interfere with the wild-type protein's function, and the phenotype would likely be recessive or haploinsufficient, not dominant-negative."
        },
        "B": {
            "is_correct": False,
            "reason": "This is incorrect. By definition, a dominant-negative mutation results in a mutant (loss-of-function) phenotype in a heterozygote, not a wild-type phenotype."
        },
        "C": {
            "is_correct": True,
            "reason": "This is the most accurate description. The mutant protein binds to (or 'aggregates' with) the wild-type protein, forming non-functional complexes. This sequestration of the functional wild-type protein into useless dimers is the classic mechanism for a dominant-negative effect, leading to an overall loss-of-function phenotype."
        },
        "D": {
            "is_correct": False,
            "reason": "This is incorrect. The term 'dominant-negative' specifically describes a type of loss-of-function mutation, not a gain-of-function mutation."
        }
    }

    # 4. Determine the logically correct answer.
    correct_choice = None
    for choice, result in analysis.items():
        if result["is_correct"]:
            correct_choice = choice
            break

    # 5. Compare the LLM's answer with the correct answer.
    if llm_answer_choice.upper() == correct_choice:
        return "Correct"
    else:
        if llm_answer_choice.upper() in analysis:
            reason_for_error = analysis[llm_answer_choice.upper()]["reason"]
            return (f"Incorrect. The provided answer '{llm_answer_choice}' is wrong. "
                    f"The correct answer is '{correct_choice}'.\n"
                    f"Reason: {reason_for_error}")
        else:
            return (f"Invalid answer choice '{llm_answer_choice}'. "
                    f"Please provide one of A, B, C, or D.")

# --- Execution of the Check ---
# The provided answer from the other LLM is "C".
llm_provided_answer = "C"
result = check_biology_question_answer(llm_provided_answer)
print(result)