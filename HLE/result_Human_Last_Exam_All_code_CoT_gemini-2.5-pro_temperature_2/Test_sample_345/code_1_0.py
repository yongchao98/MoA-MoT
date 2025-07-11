def find_museum_acquisition():
    """
    This function stores and prints the acquisition details for a specific painting.
    The information is based on established art historical records.
    """
    painting_title = "The Radionist"
    artist = "Kurt Günther"
    creation_year = 1927
    acquisition_year = 1967
    
    # Based on provenance records, the acquiring institution is known.
    acquiring_museum = "Lindenau-Museum Altenburg"

    print(f"The museum that acquired the {creation_year} tempera painting \"{painting_title}\" by {artist} in {acquisition_year} is:")
    print(acquiring_museum)

find_museum_acquisition()