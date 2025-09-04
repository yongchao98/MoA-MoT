import re

def check_correctness_of_biology_riddle(llm_response: str) -> str:
    """
    Checks the correctness of the answer to the biology riddle about protein secretion.

    This function verifies the answer by:
    1. Parsing the selected choice (A, B, C, or D) from the LLM's response.
    2. Defining the correct pathway based on a biological analysis of the riddle's clues.
    3. Comparing the selected option's pathway to the correct one and providing a detailed reason if it's incorrect.
    """
    # --- Step 1: Parse the LLM's selected answer ---
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Error: Could not find a final answer in the standard format '<<<X>>>'."
    
    selected_choice = match.group(1)

    # --- Step 2: Define the options and the correct pathway based on clues ---
    options = {
        'A': ('ribosome', 'proteasome'),
        'B': ('membrane', 'nucleus'),
        'C': ('Golgi', 'mitochondrion'),
        'D': ('cytosol', 'extracellular space')
    }

    # Clue Analysis:
    # - "ribonucleoprotein particle" (SRP) meets "nascent chain" -> Occurs in the CYTOSOL.
    # - "rough" place + "sugar" -> Protein enters the Rough ER for glycosylation.
    # - "on my way" -> Protein is transported via the secretory pathway, not degraded.
    # - Final Destination -> A primary destination for secreted proteins is the EXTRACELLULAR SPACE.
    correct_start = "cytosol"
    correct_end = "extracellular space"

    # --- Step 3: Check the selected answer and provide a reason ---
    if selected_choice == 'D':
        selected_start, selected_end = options['D']
        if selected_start == correct_start and selected_end == correct_end:
            return "Correct"
        else:
            # This case is for internal logic consistency, should not be triggered if 'D' is correct.
            return "Incorrect. The logic identified 'D' as correct, but the path defined for 'D' does not match the analysis."
    else:
        chosen_start, chosen_end = options[selected_choice]
        
        if selected_choice == 'A':
            return (f"Incorrect. The destination '{chosen_end}' is wrong. The proteasome is for protein degradation, "
                    f"which contradicts the clues of glycosylation ('sugar') and transport ('on my way') that point to the secretory pathway.")
        
        elif selected_choice == 'B':
            return (f"Incorrect. Both the start and end points are wrong. The meeting occurs in the '{correct_start}', not a generic '{chosen_start}', "
                    f"and the destination is the '{correct_end}', not the '{chosen_end}'. The process described is the secretory pathway, not nuclear import.")

        elif selected_choice == 'C':
            return (f"Incorrect. Both the start and end points are wrong. The meeting occurs in the '{correct_start}', not the '{chosen_start}' (the Golgi is a later step). "
                    f"The destination is the '{correct_end}', not the '{chosen_end}'.")
        
        else:
            # Should not be reached with the current options.
            return f"Incorrect. The selected option '{selected_choice}' does not match the biological pathway described."

# The final answer from the provided context to be checked
final_answer_from_llm = "<<<D>>>"

# Run the check
result = check_correctness_of_biology_riddle(final_answer_from_llm)
print(result)