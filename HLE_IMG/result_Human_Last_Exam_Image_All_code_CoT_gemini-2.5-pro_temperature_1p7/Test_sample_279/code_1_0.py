def classify_caries():
    """
    Analyzes the location of a dental lesion and classifies it
    using G.V. Black's classification system.
    """
    # Step 1: Define the location of the lesion from the image.
    lesion_location = "pits and fissures on the occlusal (biting) surface of a posterior tooth"
    print(f"Image analysis shows the lesion is located in the {lesion_location}.")

    # Step 2: Explain G.V. Black's classification.
    print("\nG.V. Black's Classification System:")
    print("Class I: Caries in pits and fissures on occlusal surfaces.")
    print("Class II: Caries on proximal (in-between) surfaces of posterior teeth.")
    print("Class III: Caries on proximal surfaces of anterior teeth.")
    print("Class IV: Caries on proximal surfaces of anterior teeth involving the incisal angle.")
    print("Class V: Caries on the gingival third (near the gumline) of any tooth.")
    print("Class VI: Caries on cusp tips or incisal edges.")

    # Step 3: Conclude the classification based on the location.
    # The lesion is in the pits and fissures, which matches the definition for Class I.
    final_classification_number = "I"
    
    print("\nConclusion:")
    print(f"The lesion's location matches the criteria for a Class {final_classification_number} cavity.")

# Run the classification
classify_caries()