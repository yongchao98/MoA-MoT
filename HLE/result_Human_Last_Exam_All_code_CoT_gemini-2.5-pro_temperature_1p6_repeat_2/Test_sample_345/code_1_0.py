def find_museum():
    """
    This function provides information about the acquisition of Kurt Günther's painting "The Radionist".
    The information is based on public records from the relevant museum's collection.
    """
    # Details of the artwork and its acquisition
    painting_title = "The Radionist"
    artist = "Kurt Günther"
    painting_year = 1927
    acquisition_year = 1967
    
    # The museum that acquired the painting
    museum = "Galerie Neue Meister"
    location = "Dresden, Germany"
    parent_institution = "Staatliche Kunstsammlungen Dresden (Dresden State Art Collections)"

    # Print the details in a clear statement, including the required numbers.
    print(f"The museum that acquired the {painting_year} painting '{painting_title}' by {artist} in {acquisition_year} is:")
    print(f"The {museum} in {location}.")
    print(f"It is part of the {parent_institution}.")

# Execute the function to display the answer.
find_museum()