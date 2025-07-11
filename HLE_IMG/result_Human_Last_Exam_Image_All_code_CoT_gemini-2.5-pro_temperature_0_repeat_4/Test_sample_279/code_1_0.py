def get_caries_classification(description):
    """
    Determines the G.V. Black caries classification based on a description.
    """
    classifications = {
        "Class I": "Lesions in pits and fissures on occlusal, buccal, or lingual surfaces.",
        "Class II": "Lesions on proximal surfaces of posterior teeth.",
        "Class III": "Lesions on proximal surfaces of anterior teeth.",
        "Class IV": "Lesions on proximal surfaces of anterior teeth involving the incisal angle.",
        "Class V": "Lesions on the gingival third of any tooth.",
        "Class VI": "Lesions on cusp tips or incisal edges."
    }

    # Description of the lesion seen in the image
    image_lesion_description = "pits and fissures on occlusal"

    found_class = None
    for c_class, desc in classifications.items():
        if image_lesion_description in desc.lower():
            found_class = c_class
            break

    if found_class:
        print(f"The lesion is located in the '{image_lesion_description}'.")
        print(f"According to G.V. Black's classification system, this corresponds to:")
        print(f"{found_class}: {classifications[found_class]}")
    else:
        print("Classification could not be determined from the description.")

# Run the function with the description from the image
get_caries_classification("pits and fissures on occlusal")
