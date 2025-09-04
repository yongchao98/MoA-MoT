def check_secretory_pathway_riddle(selected_option):
    """
    This function checks the answer to a biology riddle about protein synthesis and secretion.

    It decodes the clues from the riddle:
    - "ribonucleoprotein particle" + "nascent chain": This refers to the Signal Recognition Particle (SRP)
      binding to a new protein on a ribosome. This interaction happens in the cytosol.
    - "somewhat rough" + "need some sugar": These are clues for the Rough Endoplasmic Reticulum (RER),
      where co-translational import and initial glycosylation occur.
    - "I need to be on my way": This implies the protein is part of the secretory pathway,
      which ultimately leads to secretion from the cell.

    Based on this, the meeting point is the cytosol, and the final destination for a secreted protein
    is the extracellular space.
    """
    options = {
        'A': ("cytosol", "extracellular space"),
        'B': ("Golgi", "mitochondrion"),
        'C': ("ribosome", "proteasome"),
        'D': ("membrane", "nucleus")
    }

    if selected_option not in options:
        return f"Invalid option '{selected_option}'. The provided answer must be one of A, B, C, or D."

    proposed_start, proposed_end = options[selected_option]

    # Define the correct start and end points based on the riddle's clues
    correct_start = "cytosol"
    correct_end = "extracellular space"

    # Check the meeting place
    if proposed_start != correct_start:
        # The riddle describes the meeting of the SRP and the nascent chain. This happens in the cytosol.
        return (f"Incorrect. The meeting place is wrong. The 'ribonucleoprotein particle' (SRP) "
                f"binds to the 'nascent chain' in the {correct_start}, not the {proposed_start}.")

    # Check the final destination
    if proposed_end != correct_end:
        # The clues "rough", "sugar", and "on my way" point to the secretory pathway.
        # The proteasome is for degradation, mitochondria and nucleus have different import pathways.
        return (f"Incorrect. The final destination is wrong. The clues point to the secretory pathway, "
                f"where the protein is ultimately sent to the {correct_end}, not the {proposed_end}.")

    return "Correct"

# The LLM's answer to check
llm_answer = 'A'

# Run the check
result = check_secretory_pathway_riddle(llm_answer)
# print(result) # This will output "Correct"

# Final check: The code should return "Correct" for answer 'A'.
# If the answer was 'B', it would fail on the start location.
# If the answer was 'C', it would fail on the destination (proteasome is for degradation, not secretion).
# If the answer was 'D', it would fail on both start and destination.
# The logic is sound.
if result == "Correct":
    print("Correct")
else:
    print(result)