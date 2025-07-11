def solve_medical_image_query():
    """
    This function identifies the atrophied brain structure and the associated disorder from the image.
    The image shows a coronal section of a brain with marked atrophy of the caudate nucleus,
    which is a classic sign of Huntington's disease.
    """
    atrophied_structure = "caudate nucleus"
    disorder = "Huntington's disease"
    
    # The final answer is formatted as "atrophied structure, disorder".
    print(f"{atrophied_structure}, {disorder}")

solve_medical_image_query()