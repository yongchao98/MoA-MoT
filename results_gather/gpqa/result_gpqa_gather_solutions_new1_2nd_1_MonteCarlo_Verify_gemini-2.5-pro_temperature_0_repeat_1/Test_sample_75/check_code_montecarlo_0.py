def check_biology_question_answer(final_answer_text: str):
    """
    Checks the correctness of the final answer for the given biological question.

    The function encodes the key biological principles from the question as a set of rules
    and evaluates the chosen option against them.

    Principles:
    1.  **Lineage-Specific Promoter:** The red marker (mRaspberry) is only expressed AFTER an iPSC starts to differentiate. Undifferentiated iPSCs are not red.
    2.  **Apoptosis is Expected:** Apoptosis (marked by the green TUNEL stain) is a common and expected event in early embryos, especially for eliminating injected cells that fail to integrate. Therefore, a green signal is highly likely.
    3.  **Promoter Function:** A promoter controls gene expression (if/when a protein is made), not the protein's subcellular localization.
    4.  **"First Thing You Notice":** In a scientific context, this usually refers to the first *significant biological finding* that addresses the experimental question, not necessarily a trivial observation or the absolute first chronological event. The goal is to link iPSC fate (differentiation) with apoptosis.
    5.  **Key Biological Event:** A major finding in such experiments is the elimination of iPSCs that begin to differentiate inappropriately. This would be visualized as a cell that is both red (differentiating) and green (undergoing apoptosis).
    """

    # The options as presented in the final answer's prompt
    options = {
        "A": "cytoplasmic localization of the red signal",
        "B": "cell line-specific red signals label different organelles",
        "C": "there is no green signal",
        "D": "green signal colocalizes with the red signal"
    }

    # Extract the letter from the final answer format, e.g., "<<<D>>>" -> "D"
    try:
        answer_letter = final_answer_text.strip().replace("<", "").replace(">", "")
        if answer_letter not in options:
            return f"Invalid answer format. The letter '{answer_letter}' is not one of the options A, B, C, D."
    except Exception:
        return "Invalid answer format. Could not extract a letter from the input."

    # Evaluate the chosen answer based on biological principles
    if answer_letter == "A":
        return (
            "Incorrect. The answer 'cytoplasmic localization of the red signal' is a poor choice. "
            "While the mRaspberry protein is cytoplasmic, this is a static property of the reporter tool, "
            "not a dynamic experimental finding about cell fate. The question asks for a significant "
            "initial observation about the biological process, which this is not."
        )
    elif answer_letter == "B":
        return (
            "Incorrect. The answer 'cell line-specific red signals label different organelles' is biologically false. "
            "A promoter controls *if* and *when* a gene is expressed, not the subcellular destination of the resulting protein. "
            "This statement misunderstands the function of a promoter."
        )
    elif answer_letter == "C":
        return (
            "Incorrect. The answer 'there is no green signal' is biologically implausible. "
            "Apoptosis (marked by the green signal) is a fundamental process in normal embryonic development and is the "
            "expected mechanism for clearing the many injected cells that fail to integrate. A green signal is almost certain to be present."
        )
    elif answer_letter == "D":
        # This option describes the key biological event the experiment is designed to observe.
        # It shows that cells which begin to differentiate (turn red) are being eliminated by apoptosis (turn green).
        # This is the most significant and meaningful "first finding" that addresses the research question.
        return "Correct"
    else:
        # This case should not be reached if the letter is validated above, but it's good practice.
        return f"The provided answer '{answer_letter}' is not a valid option."

# The final answer provided by the LLM to be checked.
final_answer_from_llm = "<<<D>>>"

# Run the check
result = check_biology_question_answer(final_answer_from_llm)
print(result)