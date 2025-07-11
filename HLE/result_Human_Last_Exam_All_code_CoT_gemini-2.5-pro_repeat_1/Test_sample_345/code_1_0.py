def find_museum_for_painting():
    """
    This function provides information about the acquisition of the painting "The Radionist".
    The data is based on known art history records.
    """
    # Information about the artwork and its acquisition
    artist = "Kurt Günther"
    painting_title = "The Radionist"
    creation_year = 1927
    acquisition_year = 1967
    acquiring_museum = "Staatliche Kunstsammlungen Dresden"

    # Print the answer in a full sentence, including all the numbers from the query.
    print(
        f'The museum that acquired the {creation_year} tempera painting "{painting_title}" '
        f'by {artist} in {acquisition_year} is the {acquiring_museum}.'
    )

# Execute the function to print the result.
find_museum_for_painting()