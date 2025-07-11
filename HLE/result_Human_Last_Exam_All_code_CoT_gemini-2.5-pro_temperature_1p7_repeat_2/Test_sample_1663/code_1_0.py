def get_poem_info():
    """
    Provides information about a specific Spanish poem.
    """
    poet = "Rafael Alberti"
    poem_title = "Ángeles muertos (Dead Angels)"
    collection = "Sobre los ángeles (1929)"
    inspired_by_artist = "Hieronymus Bosch (El Bosco)"
    artwork_title = "The Garden of Earthly Delights"

    print(f"Poet: {poet}")
    print(f"Poem: '{poem_title}' from the collection '{collection}'")
    print(f"Inspired by: The work of painter {inspired_by_artist}")
    print(f"Specifically Directed to: The painting '{artwork_title}'")

if __name__ == "__main__":
    get_poem_info()