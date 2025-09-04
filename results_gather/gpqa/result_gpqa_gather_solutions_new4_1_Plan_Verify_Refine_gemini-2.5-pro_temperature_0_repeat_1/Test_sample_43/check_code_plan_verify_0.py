import json

def check_answer_correctness():
    """
    Checks the correctness of the final answer for the given cell biology question.

    The function encodes the biological principles behind each experimental option
    into a knowledge base. It then determines the correct answer based on this
    knowledge and compares it to the provided answer.
    """

    # Knowledge base representing the biological facts for each experimental option.
    # The question asks which experiment will NOT help.
    knowledge_base = {
        'A': {
            'description': 'Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit',
            'process_measured': 'Glucose transport',
            'location_of_process': 'Plasma Membrane',
            'is_helpful_for_mitochondria': False,
            'reasoning': 'Mitochondria do not directly take up glucose. They import pyruvate, the end product of glycolysis which occurs in the cytoplasm. The glucose uptake assay measures transport across the cell\'s plasma membrane, making it an invalid assay for isolated mitochondria.'
        },
        'B': {
            'description': "Flow cytometry after labeling with JC-1 dye",
            'process_measured': 'Mitochondrial membrane potential (ΔΨm)',
            'location_of_process': 'Mitochondria',
            'is_helpful_for_mitochondria': True,
            'reasoning': 'Mitochondrial membrane potential is a key indicator of mitochondrial health and the function of the electron transport chain. This is a standard and direct method.'
        },
        'C': {
            'description': 'Transformation of cells with recombinant luciferase and luminometer reading',
            'process_measured': 'ATP concentration',
            'location_of_process': 'Cell-wide (primarily produced in Mitochondria)',
            'is_helpful_for_mitochondria': True,
            'reasoning': 'Mitochondria are the primary producers of cellular ATP. Measuring ATP levels is a fundamental way to assess the end-product of mitochondrial function.'
        },
        'D': {
            'description': 'Confocal fluorescence microscopy after Mito-RTP staining of the cells',
            'process_measured': 'Mitochondrial-specific parameters (e.g., temperature, ROS, morphology)',
            'location_of_process': 'Mitochondria',
            'is_helpful_for_mitochondria': True,
            'reasoning': 'Mito-targeted probes allow for direct visualization and quantification of mitochondrial status, which is a valid way to assess a drug\'s effect.'
        }
    }

    # The final answer provided by the LLM analysis.
    llm_answer = 'A'

    # Determine the correct answer based on the knowledge base.
    # The correct answer is the one that is NOT helpful.
    correct_option = None
    for option, details in knowledge_base.items():
        if not details['is_helpful_for_mitochondria']:
            correct_option = option
            break

    # Check if the LLM's answer matches the correct answer.
    if llm_answer == correct_option:
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed reason.
        reason_for_correct_answer = knowledge_base[correct_option]['reasoning']
        reason_why_llm_is_wrong = f"The provided answer '{llm_answer}' is incorrect. The experiment in option '{llm_answer}' is actually helpful for studying mitochondria because {knowledge_base[llm_answer]['reasoning']}"
        
        return (f"{reason_why_llm_is_wrong}\n"
                f"The correct answer is '{correct_option}' because it describes a fundamentally flawed experiment. "
                f"Reason: {reason_for_correct_answer}")

# Execute the check and print the result.
result = check_answer_correctness()
print(result)