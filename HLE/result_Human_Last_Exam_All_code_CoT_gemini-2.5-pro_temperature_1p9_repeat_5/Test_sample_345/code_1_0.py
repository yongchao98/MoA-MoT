def find_museum_info():
    """
    This function stores and prints information about the acquisition
    of Kurt Günther's painting "The Radionist".
    The information was retrieved from art collection databases.
    """
    
    # Details of the artwork and its acquisition
    painting_year = 1927
    painting_name = "The Radionist"
    artist_name = "Kurt Günther"
    acquisition_year = 1967
    museum_name = "Staatliche Kunstsammlungen Dresden"

    # Print the answer in a full sentence
    print(
        f"The {museum_name} acquired the {painting_year} tempera painting "
        f"\"{painting_name}\" by {artist_name} in {acquisition_year}."
    )

find_museum_info()