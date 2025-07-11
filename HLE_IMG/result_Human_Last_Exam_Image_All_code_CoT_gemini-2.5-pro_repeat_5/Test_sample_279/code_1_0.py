def get_caries_classification():
    """
    Provides the caries classification for the lesion shown in the image.
    The lesion is recurrent decay around an existing restoration on the
    occlusal surface of a posterior tooth.
    """

    # G.V. Black's classification is based on the location of the caries.
    classification_system = "G.V. Black's Classification"
    lesion_location = "Pits and fissures on the occlusal surface of a posterior tooth."
    
    # Based on the location, the classification is Class I.
    class_number = 1
    classification_code = "Class I"

    print(f"Analysis based on: {classification_system}")
    print(f"Lesion Location: {lesion_location}")
    print("-" * 30)
    print(f"The caries lesion corresponds to classification code: {classification_code}")
    print(f"The number for this class is: {class_number}")
    print("-" * 30)
    print("This is because the decay, even though it is recurrent around a filling, is located in the pits and fissures of the occlusal surface.")

get_caries_classification()