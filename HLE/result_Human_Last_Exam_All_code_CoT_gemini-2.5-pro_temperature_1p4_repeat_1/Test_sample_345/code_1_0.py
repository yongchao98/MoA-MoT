def find_museum():
    """
    This function stores and prints the information about the acquisition
    of Kurt Günther's painting "The Radionist".
    """
    painting_name = "The Radionist"
    artist_name = "Kurt Günther"
    year_of_creation = 1927
    year_of_acquisition = 1967
    acquiring_museum = "Lindenau-Museum, Altenburg"

    print(
        f"The {acquiring_museum} acquired the {year_of_creation} "
        f"tempera painting \"{painting_name}\" by {artist_name} in {year_of_acquisition}."
    )

find_museum()