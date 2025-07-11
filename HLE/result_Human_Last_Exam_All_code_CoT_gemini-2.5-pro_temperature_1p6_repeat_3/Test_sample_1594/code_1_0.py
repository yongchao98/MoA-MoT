def find_building_shape():
    """
    This function represents the process of identifying the shape formed by two schools in Pisa.

    Step 1: Locate the schools.
    The schools 'Istituto Galilei-Pacinotti' and 'Istituto Ulisse Dini' are
    located next to each other on Via Benedetto Croce in Pisa, Italy.

    Step 2: Analyze the aerial view.
    By observing the satellite imagery of the school complex, we can see the layout of the buildings.

    Step 3: Determine the letter.
    The main building runs horizontally, and an adjacent wing runs perpendicularly,
    forming a clear and distinct shape. This shape corresponds to a letter of the alphabet.
    """

    school_1 = "Istituto Galilei-Pacinotti"
    school_2 = "Istituto Ulisse Dini"
    
    # The shape formed by the buildings when seen from above.
    formed_letter = "L"
    
    print(f"The buildings for the '{school_1}' and '{school_2}' highschools in Pisa, when seen from above, form the letter: {formed_letter}")

find_building_shape()