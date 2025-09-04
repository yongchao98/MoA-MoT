import re

def check_answer(llm_answer_text):
    """
    Checks the correctness of the LLM's answer to the cell biology riddle.

    The function works by:
    1.  Defining the biological facts derived from the riddle's clues.
    2.  Determining the correct start ("meeting place") and end ("destination") of the journey.
    3.  Mapping this correct journey to the provided multiple-choice options to find the correct letter.
    4.  Extracting the letter chosen by the LLM from its response.
    5.  Comparing the LLM's choice with the correct choice.
    6.  Performing a basic keyword check on the LLM's reasoning to ensure it aligns with the correct biological process.
    """
    # 1. Define biological facts from the riddle's clues
    # "ribonucleoprotein particle" + "nascent chain" + "meet" -> SRP meets protein on ribosome in cytosol
    start_location = "cytosol"
    # "rough" + "sugar" + "on my way" -> Protein enters RER, gets glycosylated, enters secretory pathway.
    # A primary destination for this pathway is secretion from the cell.
    end_location = "extracellular space"

    # 2. Define the question's options and find the correct one
    options = {
        "A": "Golgi to the mitochondrion",
        "B": "ribosome to the proteasome",
        "C": "cytosol to the extracellular space",
        "D": "membrane to the nucleus"
    }

    correct_journey_description = f"{start_location} to the {end_location}"
    correct_option_letter = None
    for letter, description in options.items():
        if description == correct_journey_description:
            correct_option_letter = letter
            break

    if correct_option_letter is None:
        # This is a sanity check for the checker itself
        return "Error in checker logic: Could not determine the correct option letter."

    # 3. Extract the final answer from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The answer format is wrong. It does not contain a final choice in the format <<<X>>>."

    llm_choice_letter = match.group(1)

    # 4. Compare the LLM's choice with the correct choice
    if llm_choice_letter != correct_option_letter:
        return (f"Incorrect: The final answer is {llm_choice_letter}, but the correct answer is {correct_option_letter}. "
                f"The riddle describes a journey from the '{start_location}' to the '{end_location}', "
                f"which corresponds to option {correct_option_letter}: '{options[correct_option_letter]}'.")

    # 5. Check if the reasoning supports the correct answer
    reasoning_text = llm_answer_text.split('<<<')[0].lower()
    required_keywords = {
        "start": ["cytosol"],
        "end": ["extracellular space", "secretion", "secreted"],
        "particle": ["srp", "signal recognition particle"],
        "location": ["rer", "rough endoplasmic reticulum"]
    }

    missing_keywords = []
    if not any(kw in reasoning_text for kw in required_keywords["start"]):
        missing_keywords.append("cytosol")
    if not any(kw in reasoning_text for kw in required_keywords["end"]):
        missing_keywords.append("extracellular space/secretion")
    if not any(kw in reasoning_text for kw in required_keywords["particle"]):
        missing_keywords.append("SRP/Signal Recognition Particle")
    if not any(kw in reasoning_text for kw in required_keywords["location"]):
        missing_keywords.append("RER/Rough Endoplasmic Reticulum")

    if missing_keywords:
        return (f"Incorrect: Although the final letter '{llm_choice_letter}' is correct, the reasoning is flawed or incomplete. "
                f"It fails to mention key concepts like: {', '.join(missing_keywords)}.")

    return "Correct"

# The provided answer from the LLM
llm_answer = """
Here is a step-by-step analysis to determine the final answer.

**1. Deconstruct the Riddle's Clues**

The riddle describes a well-known process in cell biology. Let's break it down:
*   **"the ribonucleoprotein particle says to the nascent chain as they meet"**: The "ribonucleoprotein particle" is the Signal Recognition Particle (SRP). The "nascent chain" is a new protein being synthesized on a ribosome. This meeting occurs in the **cytosol**, where protein synthesis begins.
*   **"Pause there for a minute. Let me show you in..."**: The SRP binds to the nascent chain and the ribosome, causing a temporary halt in translation. It then guides the entire complex to a specific destination.
*   **"...you really need some sugar."**: This refers to glycosylation, the process of adding sugar chains to a protein. This modification begins inside the endoplasmic reticulum (ER).
*   **"It seems somewhat rough."**: This is a direct clue for the **Rough Endoplasmic Reticulum (RER)**, which is named for its appearance due to being studded with ribosomes.
*   **"I guess this is goodbye; I need to be on my way"**: Once the protein is threaded into the RER, the SRP is released ("goodbye"). The protein is now in the secretory pathway, which, for a secreted protein, ultimately leads out of the cell.

**2. Determine the Overall Journey**

*   **Starting Point**: The meeting between the SRP and the nascent chain happens in the **cytosol**.
*   **Final Destination**: The protein enters the secretory pathway (RER -> Golgi -> vesicles). A primary destination for proteins that are glycosylated and "on their way" is secretion from the cell into the **extracellular space**.
*   **Conclusion**: The journey described is from the **cytosol to the extracellular space**.

**3. Analyze the Candidate Answers**

*   **Consensus on Rationale**: Nearly every candidate answer (16 out of 17) correctly deciphers the riddle. They all identify the key players (SRP, nascent chain), the locations (cytosol, RER), and the processes (glycosylation, secretion). There is an overwhelming consensus that the journey is from the cytosol to the extracellular space.
*   **Inconsistency in Final Letter Choice**: The final letter choices (<<<A>>>, <<<B>>>, <<<C>>>, <<<D>>>) are highly inconsistent. For example, Answers 2, 4, 6, and 11 all provide the correct rationale but conclude with <<<A>>>, while Answers 3, 5, 9, and 14 provide the same correct rationale but conclude with <<<D>>>. This indicates that the agents were likely working with different or scrambled versions of the multiple-choice options.
*   **Decision Strategy**: A simple vote count of the final letters is unreliable due to the inconsistency. The most reliable method is to trust the overwhelming consensus of the *reasoning* and match that reasoning to the specific options provided in the question.

**4. Final Conclusion**

The consensus reasoning points to a journey from the **cytosol to the extracellular space**. Matching this to the provided options:

A) Golgi to the mitochondrion
B) ribosome to the proteasome
C) **cytosol to the extracellular space**
D) membrane to the nucleus

The correct option is C.

<<<C>>>
"""

# Run the check and print the result
result = check_answer(llm_answer)
print(result)