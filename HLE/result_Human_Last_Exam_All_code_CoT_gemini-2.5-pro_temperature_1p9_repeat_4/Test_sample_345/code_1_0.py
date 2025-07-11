def find_museum():
    """
    This function stores and prints the information about the museum
    that acquired Kurt Günther's painting "The Radionist".
    """
    painting_title = "The Radionist"
    artist = "Kurt Günther"
    painting_year = 1927
    acquisition_year = 1967
    museum = "Albertinum (Staatliche Kunstsammlungen Dresden)"

    print(f'The {painting_year} tempera painting "{painting_title}" by {artist} was acquired by the {museum} in {acquisition_year}.')

if __name__ == "__main__":
    find_museum()