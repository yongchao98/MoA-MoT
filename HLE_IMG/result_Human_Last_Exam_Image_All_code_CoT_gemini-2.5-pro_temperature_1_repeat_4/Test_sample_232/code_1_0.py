def solve_medical_image_query():
    """
    This function identifies the atrophied organ and the disorder from the provided image
    and prints the result in the specified format.
    """
    # The image shows marked atrophy of the caudate nucleus.
    atrophied_structure = "Caudate nucleus"
    
    # This specific atrophy is characteristic of Huntington's disease.
    disorder = "Huntington's disease"
    
    # Print the answer in the "organ, disorder" format.
    print(f"{atrophied_structure}, {disorder}")

solve_medical_image_query()