def classify_dental_lesion():
    """
    Analyzes lesion characteristics based on an image to provide a G.V. Black classification.
    """

    # --- Image Analysis ---
    # 1. Tooth Type: Posterior tooth (molar or premolar).
    # 2. Lesion Location: Pits and fissures on the occlusal (chewing) surface.
    # 3. Additional notes: There is a pre-existing amalgam filling with recurrent (secondary)
    #    caries around its margins, as well as some primary caries in other fissures.

    lesion_location = "occlusal_pits_and_fissures"

    print("Analyzing the lesion based on its location...")
    print(f"Location identified: {lesion_location} of a posterior tooth.")

    # --- G.V. Black Classification Logic ---
    # Class I: Caries in pits or fissures on the occlusal surfaces of molars and premolars.
    # Class II: Caries on the proximal (in-between) surfaces of molars and premolars.
    # Class III: Caries on the proximal surfaces of anterior teeth (incisors/canines).
    # Class IV: Class III lesion that also involves the incisal edge.
    # Class V: Caries on the gingival (gumline) third of any tooth.
    # Class VI: Caries on cusp tips or incisal edges.

    classification = ""
    explanation = ""

    if lesion_location == "occlusal_pits_and_fissures":
        classification = "Class I"
        explanation = "This classification is used for decay located in the pits and fissures of the occlusal (chewing) surfaces of posterior teeth (molars and premolars). The image clearly shows decay in these areas, including recurrent decay around an existing restoration."
    else:
        classification = "Unknown"
        explanation = "The lesion location does not match the primary criteria for this analysis."

    print("\n--- Result ---")
    print(f"The lesion corresponds to a {classification} caries.")
    print(f"Explanation: {explanation}")

# Run the classification
classify_dental_lesion()