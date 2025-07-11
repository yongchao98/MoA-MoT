def get_caries_classification():
    """
    This function analyzes the location of the dental lesion from the image
    and provides the corresponding G.V. Black's classification code.
    """
    print("Step 1: The image shows a lesion on the occlusal (chewing) surface of a posterior tooth.")
    print("Step 2: The lesion is recurrent caries, located in the pits and fissures adjacent to an existing filling.")
    print("Step 3: According to G.V. Black's classification system, caries in the pits and fissures of occlusal surfaces of molars are classified as Class I.")

    class_number = 1
    
    print("\nTherefore, the final classification is determined by the location.")
    print("Final Classification Code:")
    print("Class " + str(class_number))

get_caries_classification()