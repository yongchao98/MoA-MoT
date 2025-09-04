import re

def check_answer(llm_answer_text):
    """
    Checks the correctness of the LLM's answer to the biology riddle.

    The riddle describes the secretory pathway:
    1.  Meeting Place: A Signal Recognition Particle (SRP, a "ribonucleoprotein particle")
        binds to a "nascent chain" on a ribosome. This occurs in the CYTOSOL.
    2.  Destination: The clues "sugar" (glycosylation in the ER), "rough" (Rough ER),
        and "on my way" (secretory pathway) point to the protein being secreted
        from the cell into the EXTRACELLULAR SPACE.
    """
    # Define the correct start and end points based on biological facts
    correct_start = "cytosol"
    correct_destination = "extracellular space"

    # Define the paths for each option
    options = {
        "A": ("cytosol", "extracellular space"),
        "B": ("ribosome", "proteasome"),
        "C": ("membrane", "nucleus"),
        "D": ("Golgi", "mitochondrion")
    }

    # Extract the final letter choice from the LLM's answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The answer is not in the required format '<<<A>>>', '<<<B>>>', etc."

    answer_key = match.group(1)
    
    if answer_key not in options:
        return f"Incorrect: The provided answer key '{answer_key}' is not a valid option."

    # Get the start and destination from the chosen option
    chosen_start, chosen_destination = options[answer_key]

    # Check if the chosen option matches the biological facts
    if chosen_start != correct_start:
        return (f"Incorrect: The chosen answer '{answer_key}' starts at '{chosen_start}'. "
                f"The riddle describes the meeting between the SRP and the nascent chain, "
                f"which occurs in the '{correct_start}'.")

    if chosen_destination != correct_destination:
        return (f"Incorrect: The chosen answer '{answer_key}' ends at '{chosen_destination}'. "
                f"The clues ('sugar', 'rough', 'on my way') indicate the protein is entering the "
                f"secretory pathway, with a final destination in the '{correct_destination}'.")

    return "Correct"

# The final answer provided by the user is <<<A>>>.
# Let's create a string that contains this answer to pass to our checker.
final_answer_from_prompt = """
The analysis and the Python script have successfully converged on a single answer by systematically applying constraints derived from the biological clues in the riddle.

1.  Initial Candidates: All four options (A, B, C, D) were considered.
2.  Constraint 1 (Meeting Place): The meeting between the "ribonucleoprotein particle" (SRP) and the "nascent chain" occurs in the cytosol. This constraint correctly eliminated options C ("membrane") and D ("Golgi") as starting points, leaving A and B.
3.  Constraint 2 (Destination): The clues "sugar" (glycosylation) and "rough" (RER) point to the secretory pathway. The ultimate destination for such a protein is often the extracellular space. This is inconsistent with the "proteasome" (for degradation), which eliminated option B.
4.  Final Candidate: The process correctly pruned the options, leaving only 'A' as the valid answer.

The logic is sound and confirmed by the checker.

<<<A>>>
"""

# Run the check
result = check_answer(final_answer_from_prompt)
print(result)