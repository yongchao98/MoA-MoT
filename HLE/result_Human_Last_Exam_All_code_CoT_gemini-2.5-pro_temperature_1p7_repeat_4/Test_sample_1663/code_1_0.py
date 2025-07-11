def find_poetry_info():
    """
    This function provides information about a specific Spanish poem.
    It identifies the poet and the work of art the poem was dedicated to.
    """
    
    # The components of the answer
    poet = "Rafael Alberti"
    artwork_title = "Carceri d'invenzione"
    artwork_translation = "Imaginary Prisons"
    artist = "Giovanni Battista Piranesi"
    poem_collection = "Sobre los ángeles (Concerning the Angels)"
    
    # Print the full, detailed answer
    print("The author of these verses is the Spanish poet:")
    print(f"Poet: {poet}")
    print("\nThe poem is part of his 1929 collection, '{}'.".format(poem_collection))
    print("\nThese verses were directed to a work of art, specifically a series of 16 engravings by the Italian artist {}.".format(artist))
    print("\nThe name of the artwork is:")
    print(f"Work of Art: {artwork_title} ({artwork_translation})")

# Execute the function to print the information
find_poetry_info()