def get_caries_classification():
    """
    Analyzes the features of a dental lesion from an image and classifies it
    according to G.V. Black's Classification system.
    """
    print("Step 1: Analyzing the lesion from the image.")
    print("- The tooth is a posterior tooth (a molar).")
    print("- There is an existing metallic (amalgam) filling.")
    print("- New decay (caries) is visible around the margin of the filling.")
    print("- The location of this decay is on the occlusal (chewing) surface, in the pits and fissures.\n")

    print("Step 2: Applying G.V. Black's Classification.")
    print("This system classifies cavities based on their location.")
    print("- Class I involves the pits and fissures on the occlusal surfaces of molars and premolars.\n")

    print("Step 3: Determining the classification code.")
    print("Since the caries in the image are on the occlusal surface of a molar, it falls into the Class I category.")
    print("This specific type of lesion, which forms around an existing restoration, is called 'recurrent' or 'secondary' caries, but its location still classifies it as Class I.")

    classification_code = "I"
    print("\nFinal Classification Code:")
    print(f"Class {classification_code}")

get_caries_classification()