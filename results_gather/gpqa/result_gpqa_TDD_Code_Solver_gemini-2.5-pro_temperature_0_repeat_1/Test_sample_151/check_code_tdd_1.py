import re

def check_llm_answer_correctness():
    """
    Checks the correctness of the provided LLM answer against the biological question.

    The question asks which protein complex would be LEAST observed in a ChIP-MS assay
    of active chromatin from G1-arrested (shmoo-forming) yeast.

    - G1 Arrest: The cell is preparing for mating, not DNA replication.
    - Active Chromatin: This implies regions of active transcription.
    - ChIP-MS: This technique will pull down proteins associated with active transcription.

    - A) Nucleosome histone complex: Found everywhere on chromatin. Abundant.
    - B) Enhancer protein complex: Directly involved in activating transcription. Abundant.
    - C) Pre-initiation complex: The core machinery of transcription. Abundant.
    - D) Pre-replication complex: Assembles in G1 for replication (in S phase), not transcription.
      While present in G1, it's not a marker of *transcriptionally active* chromatin.
      Therefore, it would be the least observed in this specific assay.

    The correct answer is D.
    """
    llm_answer = """
Excellent! The previous solution passed all tests and was confirmed as correct.

I am ready for the next question. Please provide the new problem, and I will follow the Test-Driven Development (TDD) process as before.
"""
    correct_option = "D"

    # Use regex to find a potential answer in the format <<<X>>> or a standalone letter.
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        # If the <<<X>>> format is not found, check for a standalone letter A, B, C, or D.
        # This makes the check more robust to different answer formats.
        match = re.search(r'\b([A-D])\b', llm_answer, re.IGNORECASE)

    if not match:
        return (f"The provided answer is incorrect because it does not contain a valid choice (A, B, C, or D). "
                f"The response is a generic conversational text and does not address the biological question.")

    provided_option = match.group(1).upper()

    if provided_option == correct_option:
        return "Correct"
    else:
        reasoning = (
            f"The provided answer '{provided_option}' is incorrect. The correct answer is '{correct_option}'.\n\n"
            f"Reasoning: The experiment isolates proteins from 'active chromatin' in yeast that are arrested in the G1 phase (forming a shmoo). 'Active chromatin' refers to regions undergoing active transcription.\n"
            f"- Complexes A (Nucleosome), B (Enhancer), and C (Pre-initiation) are all integral to the structure and function of active transcription. They would be highly abundant.\n"
            f"- Complex D (Pre-replication complex) assembles during G1 but its function is to prepare for DNA replication, which occurs in S phase. Since the assay targets proteins involved in transcription, not replication, the pre-replication complex would be the least observed."
        )
        return reasoning

# Execute the check and print the result
result = check_llm_answer_correctness()
print(result)