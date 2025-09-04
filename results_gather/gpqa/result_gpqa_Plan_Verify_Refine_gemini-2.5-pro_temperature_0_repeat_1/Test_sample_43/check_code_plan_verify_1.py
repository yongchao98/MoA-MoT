def check_biology_experiment_answer(answer_choice: str):
    """
    Checks the correctness of an answer for the given biology question.
    The question asks to identify the experiment that will NOT help discover a drug's mitochondrial role.
    Therefore, the correct answer is the one that does not measure a direct mitochondrial function.
    """

    # A dictionary explaining why each option is or is not a mitochondrial assay.
    experimental_options = {
        'A': {
            'is_mitochondrial_assay': True,
            'reason': "Confocal microscopy with Mito-RTP staining directly measures mitochondrial membrane potential, which is a key indicator of mitochondrial health and function."
        },
        'B': {
            'is_mitochondrial_assay': False,
            'reason': "A Glucose Uptake Assay measures the transport of glucose across the cell's plasma membrane. This is not a direct mitochondrial function. Mitochondria primarily use pyruvate (a product of glycolysis), not glucose, as fuel."
        },
        'C': {
            'is_mitochondrial_assay': True,
            'reason': "A luciferase assay measures ATP levels. Since mitochondria are the cell's primary ATP producers via oxidative phosphorylation, this is a standard method to assess overall mitochondrial function."
        },
        'D': {
            'is_mitochondrial_assay': True,
            'reason': "The specified chemical (JC-1) is a well-known dye used with flow cytometry to measure mitochondrial membrane potential. This is a direct assessment of mitochondrial status."
        }
    }

    # The correct answer is the option that is NOT a mitochondrial assay.
    correct_answer = 'B'

    # Check if the provided answer is a valid option
    if answer_choice not in experimental_options:
        return f"Invalid Answer: The choice '{answer_choice}' is not one of the possible options (A, B, C, D)."

    # Check if the provided answer is correct
    if answer_choice == correct_answer:
        return "Correct"
    else:
        # Explain why the chosen answer is incorrect
        why_choice_is_wrong = experimental_options[answer_choice]['reason']
        return (f"Incorrect. The answer '{answer_choice}' is wrong because this experiment *does* help investigate mitochondrial function. "
                f"Reason: {why_choice_is_wrong}")

# Example of how to use the checking code:
# Let's assume an LLM provided 'A' as the answer.
llm_answer = 'A'
# result = check_biology_experiment_answer(llm_answer)
# print(result)
# Expected output for 'A': "Incorrect. The answer 'A' is wrong because this experiment *does* help investigate mitochondrial function. Reason: Confocal microscopy with Mito-RTP staining directly measures mitochondrial membrane potential, which is a key indicator of mitochondrial health and function."

# Let's check the correct answer 'B'.
correct_llm_answer = 'B'
# result = check_biology_experiment_answer(correct_llm_answer)
# print(result)
# Expected output for 'B': "Correct"
