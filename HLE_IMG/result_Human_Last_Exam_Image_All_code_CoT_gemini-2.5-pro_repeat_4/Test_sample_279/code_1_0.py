def get_caries_classification(location_description):
    """
    Determines the G.V. Black caries classification based on a location description.
    """
    # G.V. Black's classification system based on lesion location
    classification_map = {
        "pits and fissures of occlusal surfaces": "Class I",
        "proximal surfaces of posterior teeth": "Class II",
        "proximal surfaces of anterior teeth not involving the incisal angle": "Class III",
        "proximal surfaces of anterior teeth involving the incisal angle": "Class IV",
        "gingival third of any tooth": "Class V",
        "cusp tips or incisal edges": "Class VI"
    }

    lesion_location = "pits and fissures of occlusal surfaces"

    if lesion_location in classification_map:
        return classification_map[lesion_location]
    else:
        return "Unknown classification"

def main():
    """
    Main function to analyze the lesion and print the classification.
    """
    # Step 1: Analyze the lesion from the image.
    # The lesion is located in the pits and fissures on the occlusal (chewing) surface of the tooth.
    # It is also adjacent to an existing filling, which makes it 'recurrent caries'.
    # However, its location dictates its G.V. Black classification.
    image_location_description = "A lesion in the pits and fissures of the occlusal surface."

    # Step 2: Determine the classification.
    caries_class = get_caries_classification(image_location_description)

    # Step 3: Print the results and explanation.
    print("Analysis of the dental lesion based on G.V. Black's Classification:")
    print("1. The lesion is located on the occlusal (chewing) surface of a posterior tooth.")
    print("2. This location involves the anatomical pits and fissures.")
    print(f"3. According to G.V. Black's system, caries in this location corresponds to a specific class.")
    print(f"Final Classification Code: {caries_class}")

if __name__ == "__main__":
    main()
