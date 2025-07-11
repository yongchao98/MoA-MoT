def find_museum():
    """
    This function stores and prints the information about the acquisition
    of the painting "The Radionist" by Kurt Günther.
    """
    # Details of the artwork and acquisition
    artist = "Kurt Günther"
    painting_title = "The Radionist"
    year_painted = 1927
    year_acquired = 1967
    
    # The acquiring institution, based on information from the
    # Dresden State Art Collections (Staatliche Kunstsammlungen Dresden, SKD).
    # The specific collection is the Galerie Neue Meister, housed in the Albertinum.
    acquiring_museum = "Galerie Neue Meister (New Masters Gallery), Albertinum, Staatliche Kunstsammlungen Dresden"

    # Print the answer in a full sentence
    print(f"In the year {year_acquired}, the {year_painted} tempera painting titled \"{painting_title}\" by {artist} was acquired by the following museum:")
    print(acquiring_museum)

if __name__ == "__main__":
    find_museum()