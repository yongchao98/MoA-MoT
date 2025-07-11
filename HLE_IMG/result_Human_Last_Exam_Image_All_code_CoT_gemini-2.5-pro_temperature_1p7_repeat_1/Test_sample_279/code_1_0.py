def get_caries_classification():
    """
    Provides the G.V. Black classification for the lesion in the image.
    
    The classification is based on the location of the tooth decay.
    - The image shows decay on the occlusal (chewing) surface of a posterior tooth.
    - Specifically, the decay is in the anatomical pits and fissures.
    - It is also secondary (or recurrent) decay, as it is adjacent to an existing filling.
    
    According to G.V. Black's classification system, caries located in pits and fissures
    on the occlusal surfaces of posterior teeth are designated as Class I.
    """
    
    class_number = 1
    
    print("Based on the G.V. Black classification system, the lesion's location corresponds to:")
    print(f"Class {class_number}")
    print("\nExplanation: This classification is used for caries occurring in the pits and fissures of the occlusal surfaces of posterior teeth (molars and premolars). The image shows secondary decay around a restoration in this precise location.")

get_caries_classification()