def find_museum():
    """
    This function provides information about the acquisition of Kurt Günther's painting.
    """
    # Information from the user's query
    painting_title = "The Radionist"
    artist_name = "Kurt Günther"
    creation_year = 1927
    acquisition_year = 1967

    # The museum that acquired the painting, based on historical records.
    acquiring_museum = "Lindenau-Museum Altenburg"

    # Print the result in a full sentence.
    print(f'The {creation_year} tempera painting "{painting_title}" by {artist_name} was acquired by the {acquiring_museum} in {acquisition_year}.')

find_museum()