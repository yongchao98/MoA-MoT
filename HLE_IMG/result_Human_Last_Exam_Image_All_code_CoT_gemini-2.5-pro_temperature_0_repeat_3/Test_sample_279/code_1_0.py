def classify_caries_lesion():
    """
    Identifies the G.V. Black classification for a caries lesion
    based on its location.
    """
    # G.V. Black's classification system
    gv_black_classes = {
        "Class I": "Lesions in pits and fissures on the occlusal surfaces of molars and premolars.",
        "Class II": "Lesions on the proximal (mesial or distal) surfaces of premolars and molars.",
        "Class III": "Lesions on the proximal surfaces of anterior teeth (incisors and canines) that do not involve the incisal angle.",
        "Class IV": "Lesions on the proximal surfaces of anterior teeth that do involve the incisal angle.",
        "Class V": "Lesions on the gingival third (cervical) of the facial or lingual surfaces of any tooth.",
        "Class VI": "Lesions on the incisal edges of anterior teeth or the cusp tips of posterior teeth."
    }

    # Description of the lesion in the provided image
    lesion_description = "The image shows a lesion with an existing restoration on the occlusal surface of a molar, located within the pits and fissures."
    target_location_feature = "pits and fissures on the occlusal surfaces"

    print("Analyzing the lesion based on G.V. Black's Classification System:")
    print("-" * 60)
    print(f"Lesion Description: {lesion_description}")
    print("-" * 60)

    # Find the matching classification
    matched_class = None
    for classification, description in gv_black_classes.items():
        if target_location_feature in description:
            matched_class = classification
            break

    if matched_class:
        print(f"The lesion's location matches the description for: {matched_class}")
        print(f"Definition: {gv_black_classes[matched_class]}")
        # The prompt asks to output the number in the final equation.
        # Since there is no equation, we will just highlight the class number.
        print("\nThe final classification code is Class I.")
        print("The Roman numeral for the class is: 1")

    else:
        print("Could not classify the lesion based on the provided information.")

classify_caries_lesion()