import re

def check_answer_correctness(question, candidate_answers, final_answer_text):
    """
    Checks the correctness of the final answer for the given biological riddle.

    Args:
        question (str): The question text.
        candidate_answers (list): A list of strings, each representing a candidate answer.
        final_answer_text (str): The final answer text provided by the LLM, including the analysis.

    Returns:
        str: "Correct" if the answer is correct, or a string explaining the error.
    """

    # Step 1: Extract the final answer choice (e.g., 'A') from the provided text.
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not find a final answer in the format <<<A>>> in the provided text."
    
    provided_answer = match.group(1)

    # Step 2: Define the correct logic based on the riddle's biological clues.
    # Clue 1: "ribonucleoprotein particle" (SRP) meets "nascent chain" (new protein on a ribosome).
    # Fact 1: This meeting occurs where protein synthesis begins, which is in the cytosol.
    correct_meeting_place = "cytosol"

    # Clue 2: "you really need some sugar" (glycosylation), "It seems somewhat rough" (Rough ER), "I need to be on my way".
    # Fact 2: This describes entry into the secretory pathway. A primary ultimate destination for this pathway is secretion from the cell.
    correct_destination = "extracellular space"

    # Step 3: Map the multiple-choice options to their biological meanings.
    options = {
        'A': ("cytosol", "extracellular space"),
        'B': ("membrane", "nucleus"),
        'C': ("ribosome", "proteasome"),
        'D': ("Golgi", "mitochondrion")
    }

    # Step 4: Determine the correct option letter based on the facts.
    correct_option_letter = None
    for option, (start, end) in options.items():
        if start == correct_meeting_place and end == correct_destination:
            correct_option_letter = option
            break
    
    if correct_option_letter is None:
        return "Error in checker logic: Could not determine the correct option based on biological facts."

    # Step 5: Compare the provided answer with the correct answer and generate a response.
    if provided_answer == correct_option_letter:
        return "Correct"
    else:
        reason = f"The provided answer '{provided_answer}' is incorrect. The correct answer is '{correct_option_letter}'.\n\n"
        reason += "Reasoning based on the riddle's constraints:\n"
        reason += f"1. **Meeting Place Constraint**: The riddle states that a 'ribonucleoprotein particle' (Signal Recognition Particle) meets a 'nascent chain' (a new protein on a ribosome). This initial interaction happens in the **cytosol**, where protein synthesis begins. \n"
        reason += f"2. **Destination Constraint**: The clues 'you really need some sugar' (glycosylation), 'It seems somewhat rough' (Rough ER), and 'I need to be on my way' indicate the protein is entering the secretory pathway. A major final destination for this pathway is secretion out of the cell into the **extracellular space**.\n\n"
        reason += f"Therefore, the correct path is from the '{correct_meeting_place}' to the '{correct_destination}', which corresponds to option '{correct_option_letter}'.\n"
        
        chosen_start, chosen_end = options.get(provided_answer, ("unknown", "unknown"))
        reason += f"The selected answer '{provided_answer}' implies a path from '{chosen_start}' to '{chosen_end}', which violates the constraints derived from the riddle."
        
        return reason

# --- You can use the code below to test the function ---

# This is the data provided in the user's prompt
question_text = '''
"Oh, I know you," the ribonucleoprotein particle says to the nascent chain as they meet. "Pause there for a minute. Let me show you in; you really need some sugar."
"It seems somewhat rough. I guess this is goodbye; I need to be on my way", the chain replies. Where did they meet, and where is the chain heading?

A) cytosol to the extracellular space
B) membrane to the nucleus
C) ribosome to the proteasome
D) Golgi to the mitochondrion
'''

candidate_answers_list = [
    # ... (all 15 candidate answers from the prompt) ...
]

final_answer_from_llm = '''
Here is a step-by-step analysis of the riddle and the provided answers.

**1. Deconstruction of the Riddle's Clues**

*   **"the ribonucleoprotein particle says to the nascent chain as they meet"**: This describes the interaction between the Signal Recognition Particle (SRP), which is a ribonucleoprotein, and a newly forming protein (a nascent polypeptide chain). This meeting occurs on a ribosome, which is located in the **cytosol** where protein synthesis begins.
*   **"Pause there for a minute. Let me show you in..."**: A key function of the SRP is to temporarily halt protein synthesis (translation) and guide the entire ribosome-nascent chain complex to the membrane of the endoplasmic reticulum.
*   **"...you really need some sugar."**: This is a direct reference to glycosylation, the process of adding sugar molecules to a protein. For proteins in the secretory pathway, this process begins inside the endoplasmic reticulum.
*   **"It seems somewhat rough."**: This is a strong clue for the **Rough Endoplasmic Reticulum (RER)**, which is named for its appearance due to being studded with ribosomes. This is the destination where the SRP docks the complex.
*   **"I guess this is goodbye; I need to be on my way"**: This signifies that the protein is entering the secretory pathway. After being synthesized into the RER, it is processed, sent to the Golgi apparatus, and then transported to its final destination. For a secreted protein, this journey ends outside the cell.

**2. Synthesis of the Pathway**

*   **Meeting Place**: The process starts in the **cytosol**.
*   **Final Destination**: The protein enters the secretory pathway (via the RER and Golgi) and is ultimately secreted from the cell into the **extracellular space**.

**3. Evaluation of the Options**

*   **A) cytosol to the extracellular space**: This perfectly matches the analysis. The meeting occurs in the cytosol, and the final destination for a secreted protein is the extracellular space.
*   **B) membrane to the nucleus**: Incorrect. The meeting is in the cytosol, not on a membrane, and the destination is not the nucleus.
*   **C) ribosome to the proteasome**: Incorrect. While the meeting involves a ribosome, the destination is not the proteasome (which is for protein degradation). The clues point to transport and modification, not destruction.
*   **D) Golgi to the mitochondrion**: Incorrect. The meeting happens long before the Golgi, and the destination is not the mitochondrion.

**4. Analysis of Candidate Answers**

A significant number of the provided answers (e.g., 2, 3, 4, 6, 7, 10, 12, 15) contain a correct, detailed analysis that logically leads to option A. However, they inexplicably conclude with a different letter (B, C, or D). This indicates a likely error in the final output step of those agents, rather than a flaw in their biological reasoning. The reasoning within the text of these answers consistently supports option A.

Conversely, the answers that selected A (e.g., 5, 8, 9, 11) have analyses that are both biologically sound and internally consistent with their final choice.

**Conclusion**

Based on a deterministic breakdown of the biological clues in the riddle, the process starts in the cytosol and the protein is headed for the extracellular space. This corresponds directly to option A. The inconsistencies in the other candidate answers reinforce this conclusion.

<<<A>>>
'''

# Run the check
result = check_answer_correctness(question_text, candidate_answers_list, final_answer_from_llm)
print(result)