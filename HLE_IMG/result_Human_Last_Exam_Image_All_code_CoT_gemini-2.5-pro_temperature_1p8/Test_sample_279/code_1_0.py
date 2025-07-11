def classify_caries_lesion():
    """
    This function explains G.V. Black's caries classification and identifies
    the classification for the lesion shown in the provided image.
    """

    classifications = {
        "Class I": "Caries affecting pits and fissures on the occlusal, buccal, or lingual surfaces of posterior teeth, or the lingual of anterior teeth.",
        "Class II": "Caries affecting the proximal (mesial or distal) surfaces of posterior teeth.",
        "Class III": "Caries affecting the proximal surfaces of anterior teeth, not involving the incisal angle.",
        "Class IV": "Caries affecting the proximal surfaces of anterior teeth that involves the incisal angle.",
        "Class V": "Caries affecting the cervical (gingival) third of the facial or lingual surfaces of any tooth.",
        "Class VI": "Caries affecting the incisal edges of anterior teeth or the cusp tips of posterior teeth."
    }

    print("G.V. Black's Caries Classification System:")
    for key, value in classifications.items():
        print(f"- {key}: {value}")

    print("\n---")
    print("Analysis of the provided image:")
    print("The image shows a carious lesion located in the pits and fissures on the occlusal (biting) surface of a molar tooth, adjacent to an existing filling.")
    print("Based on its location, this falls under the definition of Class I.")
    print("\n---")
    
    # Final answer based on the analysis
    final_classification = "Class I"
    print(f"The lesion corresponds to the caries classification code: {final_classification}")

classify_caries_lesion()