import re

def check_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer to the cell biology riddle.

    The function hardcodes the correct interpretation of the riddle and the options
    to verify the LLM's chosen answer. The riddle describes the co-translational
    translocation of a protein into the secretory pathway.

    Args:
        llm_answer_text: The full text of the LLM's response, including the final answer.

    Returns:
        "Correct" if the answer is correct, otherwise a string explaining the error.
    """
    # Step 1: Determine the correct answer based on the riddle's biological clues.
    # - "ribonucleoprotein particle" (Signal Recognition Particle or SRP) meets a "nascent chain" (new protein).
    # - This meeting occurs where protein synthesis begins on free ribosomes: the cytosol.
    correct_start = "cytosol"
    
    # - The chain is taken to a "rough" place (Rough ER) for "sugar" (glycosylation)
    #   and then goes "on my way" (entering the secretory pathway).
    # - The final destination for a secreted protein is outside the cell.
    correct_end = "extracellular space"
    
    correct_journey = (correct_start, correct_end)

    # Step 2: Define the options as presented in the question.
    options = {
        "A": ("Golgi", "mitochondrion"),
        "B": ("ribosome", "proteasome"),
        "C": ("membrane", "nucleus"),
        "D": ("cytosol", "extracellular space")
    }

    # Step 3: Extract the LLM's final choice from its response.
    # The regex looks for <<<X>>> at the very end of the string, allowing for optional whitespace.
    match = re.search(r'<<<([A-D])>>>\s*$', llm_answer_text)
    if not match:
        return "The final answer is not in the required format '<<<X>>>' at the end of the response."

    llm_choice_letter = match.group(1)

    # Step 4: Get the journey corresponding to the LLM's choice.
    llm_chosen_journey = options.get(llm_choice_letter)
    if llm_chosen_journey is None:
        # This case is unlikely with the regex used, but it's a good safeguard.
        return f"The chosen answer '{llm_choice_letter}' is not one of the valid options A, B, C, or D."

    # Step 5: Compare the LLM's chosen journey with the correct one and provide feedback.
    if llm_chosen_journey == correct_journey:
        return "Correct"
    else:
        reason = f"The chosen answer '{llm_choice_letter}' corresponds to the journey '{' to '.join(llm_chosen_journey)}'. This is incorrect. "
        if llm_chosen_journey[0] != correct_journey[0]:
            reason += f"The meeting place described in the riddle is the '{correct_start}', not '{llm_chosen_journey[0]}'. "
        if llm_chosen_journey[1] != correct_journey[1]:
            reason += f"The final destination for the protein in the secretory pathway is the '{correct_end}', not '{llm_chosen_journey[1]}'."
        return reason.strip()
