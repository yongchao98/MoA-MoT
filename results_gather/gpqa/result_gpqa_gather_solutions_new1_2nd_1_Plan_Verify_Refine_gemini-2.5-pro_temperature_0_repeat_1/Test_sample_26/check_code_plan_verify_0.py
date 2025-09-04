import re

def check_answer(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the cell biology riddle.

    Args:
        question (str): The question text containing the riddle and options.
        llm_answer_text (str): The full text of the LLM's answer, including its reasoning and final choice.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """
    # 1. Determine the correct answer based on biological facts from the riddle.
    # "ribonucleoprotein particle" (SRP) meets "nascent chain" (on a ribosome) -> This happens in the cytosol.
    correct_start = "cytosol"
    
    # "sugar" (glycosylation), "rough" (Rough ER), "on my way" (secretory pathway) -> This leads to secretion from the cell.
    correct_destination = "extracellular space"

    # 2. Parse the multiple-choice options from the question text.
    # This makes the checker robust even if the option letters change.
    options = {}
    # A simple regex to find lines like "A) Golgi to the mitochondrion"
    option_pattern = re.compile(r"([A-D])\)\s*(\w+)\s+to\s+the\s+([\w\s]+)")
    for line in question.split('\n'):
        match = option_pattern.match(line.strip())
        if match:
            letter, start, dest = match.groups()
            options[letter] = (start.strip(), dest.strip())

    if not options:
        return "Error: Could not parse the multiple-choice options from the question text."

    # 3. Find the correct letter corresponding to the correct journey.
    correct_letter = None
    for letter, (start, dest) in options.items():
        if start == correct_start and dest == correct_destination:
            correct_letter = letter
            break
    
    if correct_letter is None:
        return f"Error: The correct answer ('{correct_start}' to '{correct_destination}') was not found in the provided options."

    # 4. Extract the final answer letter from the LLM's response.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect. The final answer is not in the required format '<<<A>>>', '<<<B>>>', etc."
    
    llm_choice = match.group(1)

    # 5. Compare the LLM's choice with the correct answer.
    if llm_choice == correct_letter:
        return "Correct"
    else:
        # Provide a detailed reason for the error.
        proposed_start, proposed_dest = options.get(llm_choice, ("unknown", "unknown"))
        
        reason = f"Incorrect. The final answer is {llm_choice}, but the correct answer is {correct_letter}.\n"
        reason += f"The riddle describes a journey from the '{correct_start}' to the '{correct_destination}'.\n"
        
        if proposed_start != correct_start:
            reason += f"The chosen answer's starting point '{proposed_start}' is wrong. The meeting of the 'ribonucleoprotein particle' (SRP) and the 'nascent chain' occurs in the '{correct_start}'.\n"
        
        if proposed_dest != correct_destination:
            reason += f"The chosen answer's destination '{proposed_dest}' is wrong. The clues 'sugar' (glycosylation), 'rough' (RER), and 'on my way' (secretory pathway) indicate a final destination in the '{correct_destination}'.\n"
            
        return reason.strip()

# --- Input Data ---
# The question as provided to the LLM
question_text = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
"Oh, I know you," the ribonucleoprotein particle says to the nascent chain as they meet. "Pause there for a minute. Let me show you in; you really need some sugar."
"It seems somewhat rough. I guess this is goodbye; I need to be on my way", the chain replies. Where did they meet, and where is the chain heading?

A) Golgi to the mitochondrion
B) cytosol to the extracellular space
C) membrane to the nucleus
D) ribosome to the proteasome

Once you feel you are ready for the final answer, directly return the answer with the format <<<answer content>>> at the end of your response, e.g. <<<C>>>, <<<A>>>
"""

# The final answer provided by the LLM being checked
llm_final_answer = """
Here is a step-by-step analysis of the question, including all the careful points, to arrive at the final answer.

**Step 1: Deconstruct the riddle to identify the biological components and processes.**

*   **Careful Point 1: Identify the "characters".** The "ribonucleoprotein particle" is a classic description of the **Signal Recognition Particle (SRP)**. The "nascent chain" is a newly forming protein still being synthesized on a ribosome.
*   **Careful Point 2: Determine the meeting place.** The SRP meets the nascent chain as it emerges from a ribosome. This process begins on free ribosomes in the cell's **cytosol**. Therefore, the meeting takes place in the cytosol.

**Step 2: Analyze the dialogue to understand the journey.**

*   **Careful Point 3: Decode the journey's start.** The SRP's lines, "Pause there for a minute. Let me show you in; you really need some sugar," describe its function. It temporarily halts protein synthesis ("Pause") and guides ("show you in") the entire ribosome-protein complex to the Endoplasmic Reticulum. The mention of "sugar" is a direct reference to glycosylation, a protein modification that begins in the ER.
*   **Careful Point 4: Identify the intermediate location.** The chain's reply, "It seems somewhat rough," is an unambiguous clue for the **Rough Endoplasmic Reticulum (RER)**, which is named for its appearance due to being studded with ribosomes. This is the entry point for the secretory pathway.

**Step 3: Determine the final destination.**

*   **Careful Point 5: Understand the final path.** The chain's statement, "I guess this is goodbye; I need to be on my way," signifies that the SRP has been released and the protein is now committed to the secretory pathway. This pathway (RER → Golgi apparatus → transport vesicles) commonly leads to the protein being secreted from the cell. The final destination is therefore outside the cell, in the **extracellular space**.

**Step 4: Synthesize the complete journey and evaluate the options.**

*   **Synthesis:** The meeting occurred in the **cytosol**, and the chain is heading to the **extracellular space**.
*   **Evaluation of Options:**
    *   A) Golgi to the mitochondrion: Incorrect. The meeting is in the cytosol, not the Golgi.
    *   B) cytosol to the extracellular space: **Correct.** This perfectly matches the journey derived from the riddle's clues.
    *   C) membrane to the nucleus: Incorrect. The meeting is in the cytosol, and the destination is not the nucleus.
    *   D) ribosome to the proteasome: Incorrect. The proteasome is for protein degradation, not for a protein being glycosylated and sent "on its way".

**Analysis of Candidate Answers:**
Nearly all candidate answers correctly identify the biological pathway as being from the cytosol to the extracellular space in their detailed reasoning. The discrepancy in the final letter they choose (A, B, C, or D) is a clear indication that they were working with different lists of multiple-choice options. The overwhelming consensus of the *reasoning* supports the cytosol-to-extracellular-space pathway. Based on the options provided in this specific problem, that pathway corresponds to option B.

<<<B>>>
"""

# Run the check
result = check_answer(question_text, llm_final_answer)
print(result)