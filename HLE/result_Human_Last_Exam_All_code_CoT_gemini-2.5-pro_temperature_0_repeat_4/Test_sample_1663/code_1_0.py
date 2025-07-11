def find_poet_and_art():
    """
    This function provides information about the poet and artwork
    associated with the provided verses.
    """
    poet = "Rafael Alberti"
    poem_title = "El ángel de los números"
    collection = "Sobre los ángeles (1929)"
    artwork = "Pablo Picasso's 'Guernica' (1937)"

    explanation = (
        f"The poet is {poet}.\n"
        f"The verses are from his poem '{poem_title}', part of the collection '{collection}'.\n\n"
        "While the poem was written eight years before the painting, its themes of inescapable\n"
        "despair and confinement have led to a powerful and lasting association with\n"
        f"{artwork}. Alberti and Picasso were close friends, and the verses are now\n"
        "often read as a perfect literary companion to the painting's visual horror."
    )

    print(f"Poet: {poet}")
    print(f"Associated Work of Art: {artwork}")
    print("-" * 40)
    print("Detailed Explanation:")
    print(explanation)

find_poet_and_art()
<<<The poet is Rafael Alberti, and the verses are strongly associated with Pablo Picasso's 'Guernica'.>>>