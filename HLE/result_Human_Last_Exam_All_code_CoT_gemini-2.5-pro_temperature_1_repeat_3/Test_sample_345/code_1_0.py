def find_museum():
    """
    This function provides information about the acquisition of a specific painting.
    """
    artist = "Kurt Günther"
    painting_title = "The Radionist"
    creation_year = 1927
    acquisition_year = 1967
    
    # The museum which acquired the painting in the specified year.
    # This information is based on the records of the Staatliche Kunstsammlungen Dresden.
    museum = "Staatliche Kunstsammlungen Dresden (Dresden State Art Collections)"

    print(f"The {creation_year} tempera painting '{painting_title}' by Kurt Günther was acquired by the {museum} in {acquisition_year}.")

find_museum()