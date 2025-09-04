import re

def check_answer_correctness(question_text, llm_answer_text):
    """
    Checks the correctness of an LLM's answer for a cell biology question.

    The function encodes the scientific knowledge about the experimental methods
    and determines which one is not suitable for studying mitochondrial function.
    """

    # Scientific knowledge base about the experimental options
    # This is based on the final answer's lettering scheme.
    experiments_info = {
        'A': {
            'description': "Flow cytometry with JC-1 dye",
            'is_helpful': True,
            'reasoning': "This experiment measures mitochondrial membrane potential (ΔΨm), which is a direct and critical indicator of mitochondrial health and energy production."
        },
        'B': {
            'description': "Confocal fluorescence microscopy with Mito-RTP staining",
            'is_helpful': True,
            'reasoning': "This experiment uses a mitochondria-targeted probe (Mito-RTP) to measure a parameter related to metabolic activity (like temperature or ROS). This is a valid method to assess mitochondrial function in situ."
        },
        'C': {
            'description': "Luciferase assay for ATP",
            'is_helpful': True,
            'reasoning': "This experiment measures ATP levels. Since mitochondria are the primary producers of cellular ATP, this is a fundamental way to assess their function. Even if measuring extracellular ATP, it can be an indirect marker of cell stress/death caused by mitochondrial damage."
        },
        'D': {
            'description': "Differential centrifugation of mitochondria followed by a Glucose Uptake Assay",
            'is_helpful': False,
            'reasoning': "This experiment is fundamentally flawed. Mitochondria do not take up glucose directly from their environment. Glucose is taken up by the cell at the plasma membrane and converted to pyruvate in the cytoplasm. Pyruvate is what enters the mitochondria. Therefore, a glucose uptake assay on isolated mitochondria is a biologically nonsensical experiment."
        }
    }

    # The question asks for the experiment that will NOT help.
    # We find the key for the experiment that is not helpful.
    try:
        correct_option = next(key for key, value in experiments_info.items() if not value['is_helpful'])
    except StopIteration:
        return "Error in checker logic: No correct answer defined."

    # Extract the letter from the LLM's answer format, e.g., <<<D>>>
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Invalid answer format. Could not extract a letter from '{llm_answer_text}'."

    llm_option = match.group(1)

    # Check if the LLM's answer matches the correct option
    if llm_option == correct_option:
        return "Correct"
    else:
        # The LLM chose an option that IS helpful, but the question asked for one that is NOT.
        reason_for_incorrectness = (
            f"Incorrect. The provided answer '{llm_option}' refers to an experiment that IS helpful for discovering the mitochondrial role of a drug. "
            f"The reasoning is: {experiments_info[llm_option]['reasoning']}. "
            f"The question asks for the experiment that will NOT help. "
            f"The correct answer is '{correct_option}' because: {experiments_info[correct_option]['reasoning']}."
        )
        return reason_for_incorrectness

# The question and the final answer provided by the LLM
question = """
A eukaryotic cell evolved a mechanism to turn macromolecular building blocks into energy. The process occurs in mitochondria, which are cellular energy factories. In the series of redox reactions, the energy from food is stored between the phosphate groups and used as a universal cellular currency. The energy-laden molecules are shuttled out of the mitochondrion to serve in all cellular processes. You discovered a new anti-diabetes drug and want to investigate whether it has an effect on the mitochondria. You set up a range of experiments with your HEK293 cell line. Which of the experiments listed below will not help you discover the mitochondrial role of your drug:

A) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide
B) Confocal fluorescence microscopy after Mito-RTP staining of the cells
C) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant
D) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit
"""
final_answer = "<<<D>>>"

# Run the check
result = check_answer_correctness(question, final_answer)
print(result)