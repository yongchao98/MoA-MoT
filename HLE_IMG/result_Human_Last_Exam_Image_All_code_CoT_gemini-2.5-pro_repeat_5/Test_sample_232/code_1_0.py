def solve_medical_image_task():
    """
    This function identifies the atrophied structure and the associated disorder
    based on the provided brain scan image and prints the answer.
    """
    # The image shows marked atrophy of the caudate nucleus.
    atrophied_organ_part = "caudate nucleus"
    
    # This specific type of atrophy is a hallmark of Huntington's disease.
    disorder = "Huntington's disease"
    
    # Print the answer in the specified "organ, disorder" format.
    print(f"{atrophied_organ_part}, {disorder}")

solve_medical_image_task()