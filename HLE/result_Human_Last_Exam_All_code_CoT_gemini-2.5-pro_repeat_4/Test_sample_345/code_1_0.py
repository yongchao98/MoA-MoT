def find_museum():
    """
    This function stores and prints information about the acquisition of a specific painting.
    """
    # Known details of the artwork and its acquisition
    artist = "Kurt Günther"
    artwork_title = "The Radionist"
    creation_year = 1927
    acquisition_year = 1967
    acquiring_museum = "Lindenau-Museum Altenburg"

    # The numbers involved in the query
    num1 = creation_year
    num2 = acquisition_year

    # Construct the final answer sentence
    answer = (
        f"The {acquiring_museum} acquired the {num1} tempera painting "
        f"'{artwork_title}' by {artist} in {num2}."
    )
    print(answer)

find_museum()