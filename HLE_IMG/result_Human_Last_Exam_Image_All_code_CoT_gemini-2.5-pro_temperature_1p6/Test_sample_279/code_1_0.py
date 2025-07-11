def classify_caries_from_image_features():
    """
    This script classifies a dental lesion based on pre-identified features
    from an image, using the G.V. Black classification system.
    """
    
    # Step 1: Define the features observed in the provided image.
    # The image shows decay adjacent to a filling on the chewing surface of a molar.
    lesion_location = "Pits and fissures of the occlusal surface of a posterior tooth"
    lesion_context = "Secondary/recurrent caries adjacent to an existing restoration"

    print("Observed Lesion Features:")
    print(f"- Location: {lesion_location}")
    print(f"- Context: {lesion_context}\n")

    # Step 2: Define the G.V. Black Classification System.
    gv_black_classification = {
        "Class I": "Caries affecting pits and fissures on occlusal, buccal, and lingual surfaces of posterior teeth, and lingual of anterior teeth.",
        "Class II": "Caries affecting proximal surfaces of molars and premolars.",
        "Class III": "Caries affecting proximal surfaces of central incisors, lateral incisors, and cuspids, not involving the incisal angle.",
        "Class IV": "Caries affecting proximal including incisal angles of anterior teeth.",
        "Class V": "Caries affecting the gingival third of the facial or lingual surfaces of any tooth.",
        "Class VI": "Caries affecting cusp tips of molars, premolars, and cuspids, and incisal edges of anterior teeth."
    }
    
    print("Applying G.V. Black's Classification System...")
    
    # Step 3: Determine the classification based on the location.
    final_classification = None
    if "pits and fissures" in lesion_location.lower() and "occlusal surface" in lesion_location.lower():
        final_classification = "Class I"
        
    # Step 4: Print the result and explanation.
    if final_classification:
        print(f"\nResult: The lesion's location matches the description for {final_classification}.")
        print(f"Explanation: {gv_black_classification[final_classification]}")
        print("\nTherefore, the caries lesion corresponds to G.V. Black Class I.")
        # As requested, printing the number involved in the final classification name.
        # The Roman numeral 'I' corresponds to the number 1.
        print("Final Classification Number: 1")
    else:
        print("\nCould not determine the classification based on the provided features.")

classify_caries_from_image_features()