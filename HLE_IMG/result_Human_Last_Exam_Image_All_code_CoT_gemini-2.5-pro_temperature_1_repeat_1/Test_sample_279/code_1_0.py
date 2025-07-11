def get_caries_classification(location_description):
    """
    Determines the G.V. Black caries classification based on the lesion's location.
    """
    classifications = {
        "Class I": "Caries affecting pits and fissures on occlusal, buccal, and lingual surfaces of posterior teeth, and lingual of anterior teeth.",
        "Class II": "Caries affecting proximal surfaces of molars and premolars.",
        "Class III": "Caries affecting proximal surfaces of central incisors, lateral incisors, and cuspids not involving the incisal angle.",
        "Class IV": "Caries affecting proximal including incisal angles of anterior teeth.",
        "Class V": "Caries affecting gingival 1/3 of facial or lingual surfaces of anterior or posterior teeth.",
        "Class VI": "Caries affecting cusp tips of molars, premolars, and cuspids."
    }

    lesion_class = "Unknown"
    for c_class, desc in classifications.items():
        if "pits and fissures on occlusal" in desc and location_description == "occlusal pits and fissures":
            lesion_class = c_class
            break

    print("G.V. Black's Caries Classification System:")
    for c_class, desc in classifications.items():
        print(f"- {c_class}: {desc}")

    print("\nAnalysis:")
    print("The lesion in the image is located on the occlusal (chewing) surface, in the pits and fissures around an existing restoration.")
    print("This type of lesion is known as recurrent caries.")
    print(f"According to G.V. Black's classification, caries in this location corresponds to:")
    print(f"{lesion_class}")

# The location of the lesion in the image is the occlusal pits and fissures.
get_caries_classification("occlusal pits and fissures")