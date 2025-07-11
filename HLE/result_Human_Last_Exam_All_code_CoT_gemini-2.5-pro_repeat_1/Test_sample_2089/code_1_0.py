def find_movie_scene():
    """
    This function identifies the film and scene based on the user's query
    and prints the details.
    """
    film_title = "Midnight Cowboy"
    character_name = "Enrico 'Ratso' Rizzo"
    context = (
        "In the Oscar-winning film 'Midnight Cowboy' (1969), the character "
        f"{character_name} (played by Dustin Hoffman) is gravely ill. "
        "His friend Joe Buck (Jon Voight) takes him on a bus to Florida to fulfill his dream. "
        "After Joe buys him new clothes for the trip, a weak Ratso settles into his seat and "
        "sincerely says 'Thank you' to Joe. He dies on the bus shortly after, "
        "making this his final, fate-sealing line."
    )
    
    print(f"The film you're referring to is: {film_title}")
    print("\nScene Description:")
    print(context)

find_movie_scene()