def find_museum():
    """
    This function provides the name of the museum that acquired Kurt Günther's "The Radionist".
    """
    museum_name = "Galerie Neue Meister"
    parent_collection = "Staatliche Kunstsammlungen Dresden (Dresden State Art Collections)"
    painting = "The Radionist (Der Radionist)"
    artist = "Kurt Günther"
    year_created = 1927
    year_acquired = 1967

    print(f"The {year_created} painting '{painting}' by {artist} was acquired in {year_acquired} by the {museum_name},")
    print(f"which is part of the {parent_collection}.")

find_museum()