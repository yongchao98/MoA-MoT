def solve_medical_puzzle():
    """
    This function identifies the atrophied structure and the corresponding disorder
    based on the provided medical image.
    """

    # The image shows a coronal section of the brain, specifically the midbrain.
    # The key pathological finding is the loss of pigmentation (depigmentation)
    # in the substantia nigra, which is due to the loss of dopaminergic neurons.
    # This neuronal loss is a form of atrophy.

    # What is atrophied?
    atrophied_structure = "substantia nigra"

    # What is the disorder called?
    disorder = "Parkinson's disease"

    # The prompt requests the answer in the "organ, disorder" format.
    # We will use the specific atrophied structure for the first part of the answer.
    final_answer = f"{atrophied_structure}, {disorder}"
    
    print(final_answer)

solve_medical_puzzle()