def get_caries_classification():
    """
    Identifies the G.V. Black caries classification based on the lesion's location.
    """
    # G.V. Black's classification system based on lesion location
    gv_black_classes = {
        "Class I": "Caries affecting pits and fissures on occlusal, buccal, and lingual surfaces of posterior teeth, and lingual of anterior teeth.",
        "Class II": "Caries affecting proximal surfaces of molars and premolars.",
        "Class III": "Caries affecting proximal surfaces of central incisors, lateral incisors, and cuspids, not involving the incisal angle.",
        "Class IV": "Caries affecting proximal including incisal angles of anterior teeth.",
        "Class V": "Caries affecting the gingival 1/3 of facial or lingual surfaces of anterior or posterior teeth.",
        "Class VI": "Caries affecting cusp tips of molars, premolars, and cuspids."
    }

    # Analysis of the lesion in the image
    lesion_location_description = "The lesion in the image is on the occlusal (chewing) surface of a posterior tooth, located within the pits and fissures. There is an existing restoration with recurrent decay."

    # Determine the classification
    # Based on the location (pits and fissures of a posterior tooth), it is a Class I lesion.
    classification = "Class I"

    print("Analysis of the Dental Lesion:")
    print(f"Lesion Location: {lesion_location_description}")
    print("\nApplying G.V. Black's Classification System:")
    print(f"Definition of {classification}: {gv_black_classes[classification]}")
    print("\nConclusion:")
    print(f"The lesion corresponds to the classification code: {classification}")

get_caries_classification()