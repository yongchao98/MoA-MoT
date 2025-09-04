import re

def check_correctness_of_llm_answer(llm_answer_text: str):
    """
    Checks the correctness of the LLM's answer for the cell biology question.
    
    The question describes the co-translational targeting of a protein to the secretory pathway.
    Key clues:
    - "ribonucleoprotein particle" (SRP) meets "nascent chain" -> happens in the cytosol.
    - "sugar" and "rough" -> points to the Rough Endoplasmic Reticulum (RER).
    - "on my way" -> implies the secretory pathway (RER -> Golgi -> ...).
    - A primary destination for this pathway is secretion into the extracellular space.
    """
    
    # Extract the final letter choice from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure: Could not find a final answer in the format <<<X>>>."
    
    llm_choice = match.group(1)

    # Define the biological facts derived from the question's narrative
    # 1. Where they meet (Start)
    # The SRP (ribonucleoprotein particle) binds to the nascent chain on a ribosome.
    # This entire process begins in the cytosol.
    correct_start_location = "cytosol"

    # 2. Where the chain is heading (End)
    # The clues "rough" (Rough ER) and "sugar" (glycosylation in the ER) indicate
    # the protein has entered the secretory pathway. A major destination for proteins
    # in this pathway is secretion from the cell into the extracellular space.
    correct_end_destination = "extracellular space"

    # Define the options given in the question
    options = {
        'A': {'start': 'ribosome', 'end': 'proteasome'},
        'B': {'start': 'membrane', 'end': 'nucleus'},
        'C': {'start': 'cytosol', 'end': 'extracellular space'},
        'D': {'start': 'Golgi', 'end': 'mitochondrion'}
    }

    # Get the details of the chosen option
    chosen_option = options.get(llm_choice)
    if not chosen_option:
        return f"Failure: The chosen option '{llm_choice}' is not a valid option (A, B, C, or D)."

    # --- Verification Logic ---
    
    # Check the starting location
    # Note: While the meeting is on a ribosome, the question "Where did they meet?" refers to the
    # cellular compartment, which is the cytosol. 'cytosol' is the most accurate answer.
    if chosen_option['start'] != correct_start_location:
        return (f"Incorrect. The starting location is wrong. The question describes the meeting of the SRP and the nascent chain, "
                f"which occurs in the '{correct_start_location}'. The answer choice '{llm_choice}' incorrectly states the start is '{chosen_option['start']}'.")

    # Check the final destination
    if chosen_option['end'] != correct_end_destination:
        reason = ""
        if chosen_option['end'] == 'proteasome':
            reason = "The proteasome is for protein degradation, not a destination for a protein being glycosylated in the ER."
        elif chosen_option['end'] == 'nucleus':
            reason = "Nuclear proteins use a different import pathway and do not enter the secretory pathway."
        elif chosen_option['end'] == 'mitochondrion':
            reason = "Mitochondrial proteins use a different import pathway and do not enter the secretory pathway."
        
        return (f"Incorrect. The final destination is wrong. The clues point to the secretory pathway, with a likely destination being the '{correct_end_destination}'. "
                f"The answer choice '{llm_choice}' states the end is '{chosen_option['end']}', which is incorrect. {reason}")

    # If all checks pass, the answer is correct.
    return "Correct"