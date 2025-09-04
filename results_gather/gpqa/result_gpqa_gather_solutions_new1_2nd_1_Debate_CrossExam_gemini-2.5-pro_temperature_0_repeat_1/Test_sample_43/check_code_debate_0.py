import re

def check_answer_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the mitochondrial experiment question.

    The function uses a knowledge base to determine the biologically correct answer
    and compares it against the answer provided by the LLM.

    Args:
        llm_answer_text: A string containing the full response from the LLM,
                         which should include the final answer in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or a reason for the incorrectness.
    """

    # Step 1: Define the biological knowledge base for each experimental option.
    # The question asks which experiment will NOT help.
    knowledge_base = {
        'A': {
            'description': "Confocal fluorescence microscopy after Mito-RTP staining of the cells.",
            'is_helpful': True,
            'reasoning': "This is a direct method to visualize and quantify changes in mitochondrial activity or stress using a mitochondria-specific probe."
        },
        'B': {
            'description': "Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide (JC-1).",
            'is_helpful': True,
            'reasoning': "This is a standard and direct quantitative measure of mitochondrial membrane potential, which is a key indicator of mitochondrial health."
        },
        'C': {
            'description': "Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit.",
            'is_helpful': False,
            'reasoning': "This experiment is based on a flawed biological premise. Mitochondria do not take up glucose directly; they import pyruvate from the cytoplasm. Therefore, this experiment will not yield meaningful data."
        },
        'D': {
            'description': "Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant.",
            'is_helpful': True,
            'reasoning': "This measures ATP, the primary output of mitochondrial energy production. While measuring extracellular ATP is indirect, it's a valid indicator of cell stress/death caused by mitochondrial toxicity."
        }
    }

    # Step 2: Determine the correct answer based on the knowledge base.
    # The correct answer is the option that is NOT helpful.
    correct_option = None
    reason_for_correctness = ""
    for option, details in knowledge_base.items():
        if not details['is_helpful']:
            correct_option = option
            reason_for_correctness = details['reasoning']
            break

    if correct_option is None:
        return "Error in checker: Could not determine the correct answer from the knowledge base."

    # Step 3: Parse the provided final answer from the LLM's text.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not parse the final answer from the provided text. The answer should be in the format <<<X>>>."

    provided_option = match.group(1)

    # Step 4: Compare the provided answer with the correct answer and generate the result.
    if provided_option == correct_option:
        return "Correct"
    else:
        # The provided answer is incorrect. Explain why.
        incorrect_option_details = knowledge_base.get(provided_option)
        if incorrect_option_details and incorrect_option_details['is_helpful']:
            # The user chose a helpful experiment, but the question asked for the unhelpful one.
            return (f"Incorrect. The provided answer is '{provided_option}', but this experiment is helpful for studying mitochondria. "
                    f"The question asks for the experiment that will NOT help. "
                    f"The correct answer is '{correct_option}' because: {reason_for_correctness}")
        else:
            # This case handles if the provided answer was the wrong unhelpful one (not possible here) or other errors.
            return (f"Incorrect. The provided answer is '{provided_option}'. "
                    f"The correct answer is '{correct_option}' because: {reason_for_correctness}")
