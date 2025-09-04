import re

def check_biology_answer(question_options, final_answer_text):
    """
    Checks the correctness of a biology question answer by encoding the biological logic.

    Args:
        question_options (dict): A dictionary mapping option letters to complex names.
        final_answer_text (str): The final answer provided by the LLM, e.g., "<<<C>>>".

    Returns:
        str: "Correct" or a string explaining the error.
    """
    # Step 1 & 2: Define the cellular context and experimental goal.
    # In a yeast shmoo (G1 arrest):
    # - Primary Activity: Active Transcription
    # - Inhibited Process: DNA Replication
    # - Experimental Target: Proteins on transcriptionally active chromatin.

    # Step 3: Characterize each protein complex based on biological knowledge.
    complex_properties = {
        "enhancer protein complex": {
            "process": "Transcription",
            "status_in_shmoo": "Active",
            "reason": "Enhancers are bound by activators to drive the high levels of gene expression needed for mating. They are a hallmark of active transcription."
        },
        "nucleosome histone complex": {
            "process": "Chromatin Structure",
            "status_in_shmoo": "Structural/Ubiquitous",
            "reason": "Nucleosomes are the fundamental building blocks of all chromatin (active and inactive). They will be the most abundant protein component recovered."
        },
        "pre-replication complex": {
            "process": "DNA Replication",
            "status_in_shmoo": "Inhibited",
            "reason": "The pre-replication complex (pre-RC) prepares for DNA replication, a process that is actively blocked in G1-arrested shmoo cells. Its function is stalled."
        },
        "pre-initiation complex": {
            "process": "Transcription",
            "status_in_shmoo": "Active",
            "reason": "The pre-initiation complex (PIC) includes RNA polymerase and is essential for starting transcription. It will be abundant on the many actively transcribed mating genes."
        }
    }

    # Step 4: Identify the logically correct answer.
    # The question asks for the LEAST observed complex in an assay for ACTIVE chromatin.
    # This will be the complex associated with an inhibited process.
    least_observed_complex_name = None
    for name, properties in complex_properties.items():
        if properties["status_in_shmoo"] == "Inhibited":
            least_observed_complex_name = name
            break
    
    if not least_observed_complex_name:
        return "Error in checking logic: Could not identify the least observed complex."

    # Find the option letter corresponding to the correct complex name.
    correct_option_letter = None
    for letter, name in question_options.items():
        if name == least_observed_complex_name:
            correct_option_letter = letter
            break

    # Step 5: Verify the provided answer.
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return f"Invalid answer format. Expected '<<<X>>>' but got '{final_answer_text}'."
    
    provided_answer_letter = match.group(1)

    if provided_answer_letter == correct_option_letter:
        return "Correct"
    else:
        incorrect_complex_name = question_options.get(provided_answer_letter, "Unknown")
        incorrect_reason = complex_properties.get(incorrect_complex_name, {}).get("reason", "No reason available.")
        correct_reason = complex_properties.get(least_observed_complex_name, {}).get("reason", "No reason available.")

        return (f"Incorrect. The provided answer is {provided_answer_letter} ({incorrect_complex_name}), but this is wrong.\n"
                f"Reasoning for why {provided_answer_letter} is wrong: {incorrect_reason}\n"
                f"The correct answer is {correct_option_letter} ({least_observed_complex_name}).\n"
                f"Reasoning for why {correct_option_letter} is correct: {correct_reason}")


# --- Execution ---
# Define the options from the original question
question_options = {
    "A": "enhancer protein complex",
    "B": "nucleosome histone complex",
    "C": "pre-replication complex",
    "D": "pre-initiation complex"
}

# The final answer provided by the agent to be checked
final_answer_from_llm = "<<<C>>>"

# Run the check
result = check_biology_answer(question_options, final_answer_from_llm)
print(result)