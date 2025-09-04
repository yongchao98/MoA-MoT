import sys
import io

def check_answer_correctness():
    """
    This function checks the correctness of the provided answer to a biology riddle.
    It models the biological process described in the riddle and compares it with the given answer.
    """
    try:
        # --- Define the biological facts based on the riddle's clues ---

        # Clue 1: "the ribonucleoprotein particle says to the nascent chain as they meet."
        # This describes the Signal Recognition Particle (SRP) meeting a new protein chain on a ribosome.
        # This event occurs in the cytosol.
        meeting_location = "cytosol"

        # Clue 2 & 3: "you really need some sugar." and "It seems somewhat rough."
        # These clues point to the Rough Endoplasmic Reticulum (RER).
        # "Rough" refers to the RER's appearance.
        # "Sugar" refers to N-linked glycosylation, a modification that happens in the RER lumen.
        # This is an intermediate step in the secretory pathway.
        intermediate_location = "Rough ER"

        # Clue 4: "I guess this is goodbye; I need to be on my way"
        # This indicates the protein is being secreted from the cell.
        # The final destination for a secreted protein is outside the cell.
        final_destination = "extracellular space"

        # --- The provided answer and options ---
        # The LLM's answer is D.
        llm_answer = "D"
        options = {
            "A": ("membrane", "nucleus"),
            "B": ("ribosome", "proteasome"),
            "C": ("Golgi", "mitochondrion"),
            "D": ("cytosol", "extracellular space")
        }

        # --- Verification Logic ---
        if llm_answer not in options:
            return f"Invalid answer key '{llm_answer}'. The key must be one of {list(options.keys())}."

        chosen_start, chosen_end = options[llm_answer]

        # 1. Check the starting location ("Where did they meet?")
        # The ribosome is in the cytosol, so the meeting happens there.
        if chosen_start != meeting_location:
            return (f"Incorrect: The starting location is wrong. "
                    f"The answer claims the meeting is in the '{chosen_start}', but the interaction between a nascent chain on a ribosome and the SRP (a ribonucleoprotein particle) occurs in the '{meeting_location}'.")

        # 2. Check the final destination ("where is the chain heading?")
        # The clues about the RER and secretion point to the extracellular space.
        if chosen_end != final_destination:
            # Let's check for common misconceptions.
            if chosen_end == "proteasome":
                return (f"Incorrect: The destination is wrong. "
                        f"The answer claims the destination is the '{chosen_end}', which is for protein degradation. "
                        f"The clues about glycosylation ('sugar') and the secretory pathway ('on my way') point towards the '{final_destination}'.")
            else:
                return (f"Incorrect: The destination is wrong. "
                        f"The answer claims the destination is the '{chosen_end}'. "
                        f"However, a protein processed in the 'rough' ER and sent 'on its way' is typically secreted to the '{final_destination}'.")

        # 3. If both start and end points are correct, the answer is correct.
        return "Correct"

    except Exception as e:
        # Capture any unexpected errors during execution
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_answer_correctness()
print(result)