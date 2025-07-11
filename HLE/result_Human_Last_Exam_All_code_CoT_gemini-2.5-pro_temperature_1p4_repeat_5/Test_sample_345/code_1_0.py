def find_museum():
    """
    This function provides information about the acquisition of a specific painting.
    """
    painting_year = 1927
    painting_title = "The Radionist"
    artist = "Kurt Günther"
    acquisition_year = 1967
    acquiring_museum = "Galerie Neue Meister (part of the Staatliche Kunstsammlungen Dresden)"

    print(f"The {painting_year} tempera painting \"{painting_title}\" by {artist} was acquired in {acquisition_year} by the {acquiring_museum}.")

if __name__ == "__main__":
    find_museum()