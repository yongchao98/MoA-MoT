def solve_poetry_riddle():
    """
    This function identifies the poet and artwork for the given verses
    and prints the answer, including relevant years in an 'equation' format.
    """
    poet = "Rafael Alberti"
    artwork_name = "Guernica"
    artist = "Pablo Picasso"
    poem_collection = "A la pintura"
    painting_year = 1937
    collection_publication_year = 1945

    print(f"The poet who wrote these verses is {poet}.")
    print(f"They were directed to the work of art: '{artwork_name}' by {artist}.")
    print(f"\nThe verses are from his poetry collection '{poem_collection}'.")

    # The final "equation" displaying the numbers
    print("\n--- Contextual Equation ---")
    print(f"Painting Year ({painting_year}) + Poem Publication Year ({collection_publication_year})")


solve_poetry_riddle()