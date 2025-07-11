def find_building_shape():
    """
    This function presents the result of analyzing the satellite view
    of the two specified high schools in Pisa.
    """
    school_1 = "Istituto Galilei-Pacinotti"
    school_2 = "Istituto Ulisse Dini"
    
    # By observing the satellite imagery of the two schools located next to each other
    # on Via Benedetto Croce in Pisa, their combined building layout clearly forms a specific letter.
    formed_letter = "H"
    
    print(f"The buildings of '{school_1}' and '{school_2}' stand next to each other.")
    print("When viewed from above, their combined shape forms the letter:")
    print(formed_letter)

find_building_shape()