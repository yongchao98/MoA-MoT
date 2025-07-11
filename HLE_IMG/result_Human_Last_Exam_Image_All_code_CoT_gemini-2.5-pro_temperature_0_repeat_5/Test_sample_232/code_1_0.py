def solve_medical_case():
    """
    This function identifies the atrophied structure and the disorder based on the provided image.
    The image shows a coronal section of a brain.
    There is marked atrophy of the head of the caudate nucleus, which normally bulges into the lateral ventricles.
    This specific atrophy is a hallmark of Huntington's disease.
    """
    
    # The specific part of the organ that has atrophied
    atrophied_structure = "caudate nucleus"
    
    # The name of the disorder
    disorder_name = "Huntington's disease"
    
    # Print the answer in the specified format
    print(f"{atrophied_structure}, {disorder_name}")

solve_medical_case()