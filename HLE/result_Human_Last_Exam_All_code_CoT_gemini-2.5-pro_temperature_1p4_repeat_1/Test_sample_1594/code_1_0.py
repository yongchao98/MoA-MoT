def find_building_shape():
    """
    This function explains the shape formed by the two high schools.
    
    The problem requires identifying the letter shape formed by the buildings of:
    1. Istituto Galilei-Pacinotti
    2. Istituto Ulisse Dini
    
    By looking at these schools on a satellite map, we can see their layout.
    The Liceo Dini building forms a long vertical line. The Galilei-Pacinotti
    building forms a horizontal line across the top of the Dini building.
    
    When combined, they form a specific letter.
    """
    
    # The letter formed by the buildings
    formed_letter = "T"
    
    print(f"The buildings of Istituto Galilei-Pacinotti and Istituto Ulisse Dini in Pisa, when seen from above, form the letter: {formed_letter}")

if __name__ == "__main__":
    find_building_shape()