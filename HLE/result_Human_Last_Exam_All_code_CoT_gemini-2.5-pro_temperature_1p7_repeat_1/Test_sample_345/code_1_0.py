def find_museum():
    """
    This function provides the name of the museum that acquired Kurt Günther's painting "The Radionist".
    """
    painting_name = "The Radionist"
    artist_name = "Kurt Günther"
    creation_year = 1927
    acquisition_year = 1967
    
    # Based on art historical records, the museum is the Lindenau-Museum Altenburg.
    museum_name = "Lindenau-Museum Altenburg"
    
    print(f"The museum that acquired the {creation_year} painting '{painting_name}' by {artist_name} in {acquisition_year} is:")
    print(museum_name)

find_museum()