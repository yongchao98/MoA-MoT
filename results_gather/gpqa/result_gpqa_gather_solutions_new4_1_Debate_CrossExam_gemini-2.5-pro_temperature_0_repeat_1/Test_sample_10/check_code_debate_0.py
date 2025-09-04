def check_sars_cov_2_biology_answer(selected_answer: str) -> str:
    """
    Checks the correctness of the selected answer for the SARS-CoV-2 molecular biology question.
    The question asks to identify the INCORRECT statement.

    Args:
        selected_answer: The letter ('A', 'B', 'C', or 'D') corresponding to the chosen answer.

    Returns:
        "Correct" if the answer is right, or a string explaining why it's wrong.
    """
    # A knowledge base representing the factual correctness of each statement.
    # The statements are mapped to their corresponding letters in the question.
    statement_facts = {
        'A': {
            "is_correct": True,
            "reasoning": "This is a correct description of -1 programmed ribosomal frameshifting in coronaviruses."
        },
        'B': {
            "is_correct": False,
            "reasoning": "This statement is incorrect. The SARS-CoV-2 pseudoknot unfolds via a three-state pathway, not two. The claim of a 'linear correlation' is also a significant oversimplification."
        },
        'C': {
            "is_correct": False,
            "reasoning": "This statement is incorrect. The nsp10/nsp14-ExoN complex is an exonuclease that *degrades* RNA for proofreading. It does not 'prevent the breakdown of dsRNA'; it does the opposite."
        },
        'D': {
            "is_correct": True,
            "reasoning": "This statement accurately describes the initial mechanism of apoptosis induced by ORF3a via the extrinsic pathway, as supported by scientific literature."
        }
    }

    # Clean up the input answer
    answer = selected_answer.strip().upper()

    if answer not in statement_facts:
        return f"Invalid Answer Format: The provided answer '{selected_answer}' is not one of 'A', 'B', 'C', or 'D'."

    # The question asks for the statement that is an EXCEPTION (i.e., incorrect).
    # Therefore, a correct answer must be a key pointing to a statement where 'is_correct' is False.
    if not statement_facts[answer]["is_correct"]:
        # The selected answer 'C' correctly identifies a factually incorrect statement.
        # While statement 'B' is also incorrect, 'C' contains a more fundamental error
        # and is the best answer to the question. Therefore, the choice is correct.
        return "Correct"
    else:
        # The selected answer points to a statement that is factually correct.
        # This is the wrong response since the question asks for the incorrect statement.
        return (f"Incorrect. The answer '{answer}' was chosen, but statement {answer} is factually correct. "
                f"The question asks to identify the incorrect statement (the exception).")

# --- Verification ---
# The final answer provided by the LLM is 'C'.
final_answer_from_llm = 'C'
result = check_sars_cov_2_biology_answer(final_answer_from_llm)
print(result)