def get_caries_classification(image_description):
    """
    Determines the G.V. Black caries classification based on a description.

    Args:
        image_description (str): A description of the lesion's location.

    Returns:
        str: The classification and its definition.
    """
    classifications = {
        "Class I": "Caries affecting pits and fissures on the occlusal, buccal, and lingual surfaces of posterior teeth, and the lingual of anterior teeth.",
        "Class II": "Caries affecting the proximal surfaces of premolars and molars.",
        "Class III": "Caries affecting the proximal surfaces of incisors and canines that do not involve the incisal angle.",
        "Class IV": "Caries affecting the proximal surfaces of incisors and canines that do involve the incisal angle.",
        "Class V": "Caries affecting the cervical (gingival) third of the facial or lingual surfaces of any tooth.",
        "Class VI": "Caries affecting the incisal edges of anterior teeth and the cusp tips of posterior teeth."
    }

    lesion_class = "Unknown"
    # The lesion in the image is on the occlusal surface in the pits and fissures.
    if "occlusal" in image_description and "pits and fissures" in image_description:
        lesion_class = "Class I"

    if lesion_class != "Unknown":
        print(f"Image Analysis: The lesion is located in the {image_description}.")
        print("-" * 20)
        print(f"According to G.V. Black's classification system, this corresponds to: {lesion_class}")
        print(f"Definition: {classifications[lesion_class]}")
    else:
        print("Could not classify the lesion based on the provided description.")

# Description of the lesion from the provided image.
lesion_location = "pits and fissures of the occlusal surface"
get_caries_classification(lesion_location)