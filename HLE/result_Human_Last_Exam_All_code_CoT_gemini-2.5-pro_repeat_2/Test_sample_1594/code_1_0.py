def find_building_shape():
    """
    This function identifies the letter formed by the adjacent buildings of
    two high schools in Pisa.
    """
    school_1 = "Istituto Galilei-Pacinotti"
    school_2 = "Istituto Ulisse Dini"
    
    # By observing the satellite view of the two schools, which are two
    # parallel buildings connected by a central structure, we can identify the letter.
    formed_letter = "H"
    
    print(f"The buildings for {school_1} and {school_2} stand next to each other.")
    print(f"When viewed from above, their combined shape forms the letter: {formed_letter}")

find_building_shape()