def find_museum():
    """
    This function stores and prints the information about the museum
    that acquired Kurt Günther's painting "The Radionist".
    """
    artist = "Kurt Günther"
    painting_name = "The Radionist"
    creation_year = 1927
    acquisition_year = 1967
    museum = "Staatliche Kunstsammlungen Dresden (Dresden State Art Collections)"

    print(f"The museum that acquired the {creation_year} painting '{painting_name}' by {artist} in {acquisition_year} is:")
    print(museum)

find_museum()