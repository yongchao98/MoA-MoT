def find_museum():
    """
    This function stores and prints the information about the museum
    that acquired Kurt Günther's painting "The Radionist".
    """
    painting_name = "The Radionist"
    artist_name = "Kurt Günther"
    creation_year = 1927
    acquisition_year = 1967
    
    # According to online art collection databases, the Lindenau-Museum Altenburg
    # acquired the painting in 1967.
    acquiring_museum = "Lindenau-Museum Altenburg"
    
    print(f"The {acquiring_museum} acquired the {creation_year} tempera painting \"{painting_name}\" by {artist_name} in {acquisition_year}.")

find_museum()