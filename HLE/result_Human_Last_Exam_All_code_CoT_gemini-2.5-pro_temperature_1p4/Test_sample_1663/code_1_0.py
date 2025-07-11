def solve_poetry_query():
    """
    This function identifies the poet and artwork related to the provided verses
    and prints the information in a structured way.
    """
    poet = "Rafael Alberti"
    poem = "Los ángeles muertos (The Dead Angels)"
    collection = "Sobre los ángeles (Concerning the Angels)"
    artist = "Francisco de Goya"
    artwork = "Los Caprichos (The Caprices)"

    print("The provided verses were written by the Spanish poet Rafael Alberti.")
    print("-" * 20)
    print(f"Poet: {poet}")
    print(f"Poem: '{poem}' from the collection '{collection}' (1929).")
    print("\nThe poem was not directed at a single painting, but is widely considered to be inspired by a famous series of etchings:")
    print(f"Artist: {artist}")
    print(f"Artwork: The series of 80 etchings titled '{artwork}'")

solve_poetry_query()