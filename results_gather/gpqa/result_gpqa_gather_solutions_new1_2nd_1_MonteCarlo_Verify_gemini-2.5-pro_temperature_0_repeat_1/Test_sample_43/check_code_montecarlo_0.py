def check_correctness():
    """
    This function checks the correctness of the answer to the biological question.
    It encodes the biological validity of each experimental option and compares it
    against the provided answer.
    """
    # The question asks to identify the experiment that will NOT help discover the mitochondrial role of a drug.
    # We encode the biological facts for each option.
    experimental_validity = {
        'A': {
            'description': "Flow cytometry with JC-1 dye to measure mitochondrial membrane potential.",
            'is_helpful': True,
            'reason': "Measuring mitochondrial membrane potential is a direct and standard method to assess mitochondrial health and function."
        },
        'B': {
            'description': "Confocal microscopy with a mitochondria-specific stain (Mito-RTP).",
            'is_helpful': True,
            'reason': "Using a mitochondria-specific probe to measure parameters like ROS, temperature, or morphology is a direct way to assess mitochondrial activity and stress."
        },
        'C': {
            'description': "Glucose Uptake Colorimetric Assay on isolated mitochondria.",
            'is_helpful': False,
            'reason': "This experiment is fundamentally flawed. Mitochondria do not take up glucose; they take up pyruvate, the end-product of glycolysis which occurs in the cytoplasm. Therefore, this experiment will not yield useful information."
        },
        'D': {
            'description': "Luciferase assay to measure ATP in the supernatant.",
            'is_helpful': True,
            'reason': "Measuring ATP, the primary output of mitochondria, is a valid way to assess cellular energy status. Measuring it in the supernatant can indicate cell death caused by mitochondrial toxicity."
        }
    }

    # The final answer provided by the LLM analysis.
    provided_answer = 'C'

    # The question asks for the option that is NOT helpful.
    # We need to find the key in our dictionary where 'is_helpful' is False.
    correct_key = None
    for key, value in experimental_validity.items():
        if not value['is_helpful']:
            correct_key = key
            break

    if provided_answer == correct_key:
        return "Correct"
    else:
        # If the provided answer is wrong, explain why.
        if provided_answer in experimental_validity:
            # This case handles when the user selects a valid experiment, but the question asks for the invalid one.
            if experimental_validity[provided_answer]['is_helpful']:
                return (f"Incorrect. The answer '{provided_answer}' corresponds to an experiment that IS helpful for studying mitochondria. "
                        f"The question asks for the experiment that will NOT help. "
                        f"The correct answer is '{correct_key}' because: {experimental_validity[correct_key]['reason']}")
            else:
                # This case is unlikely but handles logical inconsistencies.
                return f"Incorrect. The logic indicates '{correct_key}' is the correct answer, but '{provided_answer}' was given."
        else:
            return f"Incorrect. The provided answer '{provided_answer}' is not one of the valid options (A, B, C, D)."

# Run the check and print the result.
result = check_correctness()
print(result)