import re

def check_answer_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the given biology riddle.

    The riddle describes the co-translational translocation of a protein into the ER,
    which is the start of the secretory pathway.

    1.  **Meeting Place**: A "ribonucleoprotein particle" (Signal Recognition Particle, SRP)
        meets a "nascent chain" (new protein on a ribosome). This happens in the **cytosol**.
    2.  **Destination**: The protein is guided to the "rough" ER, gets "sugar" (glycosylated),
        and is sent "on its way". This describes the secretory pathway. The ultimate destination
        for a secreted protein is outside the cell, in the **extracellular space**.

    Therefore, the correct journey is from the cytosol to the extracellular space.
    """
    # Define the correct start and end locations based on the riddle's clues.
    correct_start = "cytosol"
    correct_end = "extracellular space"

    # Map the options to their respective start and end locations.
    options = {
        "A": ("ribosome", "proteasome"),
        "B": ("Golgi", "mitochondrion"),
        "C": ("membrane", "nucleus"),
        "D": ("cytosol", "extracellular space")
    }

    # Determine the correct option letter.
    correct_option_letter = None
    for option, (start, end) in options.items():
        # Option D perfectly matches the derived start and end points.
        if start == correct_start and end == correct_end:
            correct_option_letter = option
            break

    # Extract the chosen option from the LLM's answer.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect format: The final answer is not enclosed in '<<<>>>' or is missing."

    llm_choice = match.group(1)

    # Check if the LLM's choice is the correct one.
    if llm_choice == correct_option_letter:
        return "Correct"
    else:
        # If incorrect, provide a detailed reason.
        chosen_start, chosen_end = options[llm_choice]
        reason = f"The answer '{llm_choice}' is incorrect.\n"
        
        # Analyze the starting point.
        if chosen_start != correct_start:
            reason += f"- Constraint not satisfied: The meeting place is wrong. The riddle describes the Signal Recognition Particle meeting the nascent chain on a ribosome, a process that occurs in the '{correct_start}', not the '{chosen_start}'.\n"
        
        # Analyze the destination.
        if chosen_end != correct_end:
            reason += f"- Constraint not satisfied: The destination is wrong. The clues ('rough' ER, 'sugar' addition, 'on my way') point to the secretory pathway. A primary destination for this pathway is the '{correct_end}', not the '{chosen_end}'. The '{chosen_end}' is incorrect because "
            if chosen_end == "proteasome":
                reason += "the proteasome is for protein degradation, not secretion."
            elif chosen_end == "mitochondrion":
                reason += "mitochondria have their own distinct protein import system."
            elif chosen_end == "nucleus":
                reason += "nuclear import is a separate pathway and does not involve the rough ER."
            else:
                reason += "it does not fit the secretory pathway described."

        return reason.strip()

# The provided answer from the prompt to be checked.
llm_answer = """
Here is a step-by-step analysis of the question:

1.  **Identify the characters and setting:**
    *   **"the ribonucleoprotein particle"**: This is the Signal Recognition Particle (SRP). It is a complex made of RNA and protein, hence "ribonucleoprotein". Its job is to recognize and target specific proteins for secretion or insertion into membranes.
    *   **"the nascent chain"**: This is a newly forming (nascent) polypeptide chain that is being synthesized by a ribosome.
    *   **"as they meet"**: The SRP meets the nascent chain while it is still attached to the ribosome. This process of protein synthesis (translation) begins on free ribosomes in the **cytosol**. Therefore, the meeting place is the cytosol.

2.  **Analyze the dialogue for clues about the process:**
    *   **"Pause there for a minute."**: When the SRP binds to the nascent chain's signal sequence, it causes a temporary halt in translation. This is a key feature of the SRP-mediated pathway.
    *   **"Let me show you in; you really need some sugar."**: The SRP guides the entire ribosome-nascent chain complex to a specific location. The phrase "you really need some sugar" refers to glycosylation, the process of adding sugar molecules to the protein. This modification begins in the Endoplasmic Reticulum (ER).
    *   **"It seems somewhat rough."**: This is a direct and strong clue for the **Rough Endoplasmic Reticulum (RER)**. The RER gets its "rough" appearance from the many ribosomes studded on its surface. The SRP docks the ribosome complex onto the RER membrane.
    *   **"I guess this is goodbye; I need to be on my way"**: Once the nascent chain is successfully handed off to a channel in the RER membrane, the SRP is released ("goodbye") and recycled. The protein is now "on its way" through the secretory pathway.

3.  **Determine the final destination:**
    *   The protein has entered the RER, the first step of the secretory pathway. From the RER, it will travel to the Golgi apparatus for further processing and packaging. For a protein that is to be secreted from the cell, it will be packaged into vesicles that fuse with the plasma membrane, releasing the protein outside the cell.
    *   Therefore, the ultimate destination for this protein is the **extracellular space**.

4.  **Synthesize the answer:**
    *   **Where did they meet?** In the **cytosol**.
    *   **Where is the chain heading?** To the **extracellular space**.

5.  **Match with the given options:**
    *   A) ribosome to the proteasome: Incorrect. The proteasome is for degradation.
    *   B) Golgi to the mitochondrion: Incorrect. The meeting is in the cytosol, not the Golgi.
    *   C) membrane to the nucleus: Incorrect. The meeting is in the cytosol, and the pathway does not lead to the nucleus.
    *   D) cytosol to the extracellular space: Correct. The meeting is in the cytosol, and the final destination is the extracellular space.

<<<D>>>
"""

# Execute the check
result = check_answer_correctness(llm_answer)
print(result)