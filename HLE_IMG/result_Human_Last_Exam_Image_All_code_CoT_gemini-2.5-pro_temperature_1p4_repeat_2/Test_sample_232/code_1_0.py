def solve_disease_riddle():
    """
    This function identifies the atrophied organ part and the disorder from the image provided.
    The image shows a coronal section of a brain with significant atrophy of the caudate nucleus,
    which is characteristic of Huntington's disease.
    """
    atrophied_structure = "Caudate nucleus"
    disorder = "Huntington's disease"
    
    # The final answer format is "organ, disorder"
    final_answer = f"{atrophied_structure}, {disorder}"
    
    print(final_answer)

solve_disease_riddle()