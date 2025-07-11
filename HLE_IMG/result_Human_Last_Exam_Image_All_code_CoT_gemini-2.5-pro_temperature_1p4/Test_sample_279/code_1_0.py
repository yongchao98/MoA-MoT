def find_caries_classification():
    """
    This script determines the G.V. Black caries classification for the lesion
    shown in the image.
    """

    # Step 1: Define G.V. Black's classifications for context.
    print("Step 1: G.V. Black's Classification System")
    print("---------------------------------------------")
    print("Class I:   Caries in pits or fissures on occlusal surfaces of molars/premolars.")
    print("Class II:  Caries on proximal (in-between) surfaces of molars/premolars.")
    print("Class III: Caries on proximal surfaces of incisors/canines (front teeth).")
    print("Class IV:  Caries on proximal surfaces of incisors/canines, involving the incisal edge.")
    print("Class V:   Caries on the gingival (gumline) third of any tooth.")
    print("---------------------------------------------")

    # Step 2: Analyze the provided image.
    print("\nStep 2: Image Analysis")
    print("---------------------------------------------")
    print("The image displays the top (occlusal) surface of a molar tooth.")
    print("A large filling is present in the center, and there is dark staining and decay in the surrounding grooves (pits and fissures).")
    print("---------------------------------------------")

    # Step 3: Match the analysis to the classification.
    print("\nStep 3: Conclusion")
    print("---------------------------------------------")
    lesion_location = "Pits and fissures on the occlusal surface"
    classification_code = "I" # Roman numeral for 1

    print(f"The lesion is located in the '{lesion_location}'.")
    print(f"According to G.V. Black's system, this corresponds to a Class {classification_code} lesion.")
    print("---------------------------------------------")


# Execute the function to get the answer.
find_caries_classification()