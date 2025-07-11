def solve_medical_image_query():
    """
    This function identifies the atrophied organ part and the disorder shown in the image.
    The image shows a coronal section of a brain.
    There is visible atrophy of the caudate nucleus, which is a classic sign of Huntington's disease.
    The function will print the answer in the format "organ, disorder".
    """
    # The specific part of the organ that has atrophied.
    atrophied_part = "Caudate nucleus"

    # The name of the disorder.
    disorder = "Huntington's disease"

    # Construct the final answer string.
    answer = f"{atrophied_part}, {disorder}"

    # Print the final answer.
    print(answer)

solve_medical_image_query()