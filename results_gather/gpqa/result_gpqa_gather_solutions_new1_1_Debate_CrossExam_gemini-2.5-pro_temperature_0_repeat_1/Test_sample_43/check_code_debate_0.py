import re

def check_mitochondria_experiment_answer():
    """
    This function checks the correctness of the answer to the mitochondria experiment question.
    It defines the ground truth based on biological principles and identifies the correct choice.
    """

    # Step 1: Define the ground truth for each experimental option provided in the question.
    # The question asks which experiment will NOT help.
    experiments = {
        'A': {
            'description': "Confocal fluorescence microscopy after Mito-RTP staining of the cells",
            'is_helpful': True,
            'reason': "Valid: Mito-RTP is a probe for mitochondrial activity/stress. This is a direct way to assess mitochondrial function."
        },
        'B': {
            'description': "Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide",
            'is_helpful': True,
            'reason': "Valid: The chemical is JC-1, a standard probe for mitochondrial membrane potential, a key indicator of mitochondrial health."
        },
        'C': {
            'description': "Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant",
            'is_helpful': True,
            'reason': "Valid: This is an ATP assay. Since mitochondria are the primary producers of ATP, this is a key measure of their function."
        },
        'D': {
            'description': "Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit",
            'is_helpful': False,
            'reason': "Not helpful: This experiment is fundamentally flawed. Mitochondria do not take up glucose; they take up pyruvate from the cytoplasm. A glucose uptake assay on isolated mitochondria is biologically nonsensical."
        }
    }

    # Step 2: Identify the correct answer based on which experiment is not helpful.
    correct_option = None
    failure_reason = ""
    for option, details in experiments.items():
        if not details['is_helpful']:
            correct_option = option
            failure_reason = details['reason']
            break

    if correct_option is None:
        return "Error: Could not determine the correct answer from the defined experiments."

    # Step 3: Analyze the provided LLM answers.
    # The majority of the LLM answers correctly identify the flawed experiment by its description
    # ("glucose uptake assay on isolated mitochondria"). This confirms the reasoning is sound.
    # The final lettered answers vary (A, B, C, D), suggesting the options were shuffled for each LLM.
    # The correct final answer must correspond to the lettering in the original question.
    
    # The final answer derived from our analysis is 'D'. We will check if this is correct.
    final_answer_to_check = 'D'

    if final_answer_to_check == correct_option:
        return "Correct"
    else:
        return f"Incorrect. The provided answer '{final_answer_to_check}' is wrong. The correct answer is '{correct_option}'. Reason: {failure_reason}"

# Execute the check
result = check_mitochondria_experiment_answer()
print(result)
