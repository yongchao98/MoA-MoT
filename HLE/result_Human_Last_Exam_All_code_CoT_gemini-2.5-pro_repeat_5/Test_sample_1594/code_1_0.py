def find_building_shape():
    """
    This function provides the answer to a knowledge-based question
    about the shape of school buildings in Pisa.
    """
    school_1 = "Istituto Galilei-Pacinotti"
    school_2 = "Istituto Ulisse Dini"
    
    # Based on satellite imagery, the two adjacent school buildings
    # form the shape of a specific letter.
    formed_letter = "H"
    
    print(f"The buildings of the '{school_1}' and '{school_2}' highschools in Pisa, when seen from above, form the letter:")
    print(formed_letter)

find_building_shape()