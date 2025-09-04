def check_answer_correctness(llm_answer: str):
    """
    Checks the correctness of the given answer for designing a retrovirus diagnostic kit.

    Args:
        llm_answer: The answer provided by the LLM (e.g., 'A', 'B', 'C', 'D').

    Returns:
        A string indicating "Correct" or the reason for the error.
    """
    correct_answer = 'B'
    
    # Define the reasons why other options are incorrect based on scientific principles.
    reasons_for_incorrectness = {
        'A': "This answer is incorrect. While ELISA is a valid diagnostic method, targeting IgG antibodies is not suitable for *quick* diagnosis in a new outbreak. IgG is an indicator of a later-stage or past infection, not an early, active one.",
        'C': "This answer is incorrect. Identifying a virus based on patient symptoms is not a reliable or specific method for designing a molecular diagnostic kit. Different viruses can cause similar symptoms, and a molecular test requires precise genetic information.",
        'D': "This answer is incorrect. A retrovirus has an RNA genome, not a DNA genome. Therefore, the first step cannot be DNA sequencing of the virus. The process must involve converting the viral RNA to cDNA using reverse transcriptase."
    }

    if llm_answer.strip().upper() == correct_answer:
        # The answer is B. Let's double-check if it satisfies all constraints.
        # 1. Handles RNA genome? Yes, 'cDNA sequencing' implies reverse transcription from RNA.
        # 2. Is it quick? Yes, 'real time PCR' is a gold standard for rapid detection.
        # 3. Is it specific? Yes, sequencing provides the specific genetic information needed.
        # All constraints are met.
        return "Correct"
    elif llm_answer.strip().upper() in reasons_for_incorrectness:
        # The answer is one of the incorrect options.
        return reasons_for_incorrectness[llm_answer.strip().upper()]
    else:
        return f"Invalid answer choice provided: '{llm_answer}'. The answer must be A, B, C, or D."

# The answer provided by the LLM
llm_provided_answer = "B"

# Run the check
result = check_answer_correctness(llm_provided_answer)
print(result)