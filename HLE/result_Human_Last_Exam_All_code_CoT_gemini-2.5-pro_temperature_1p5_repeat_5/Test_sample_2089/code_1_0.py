def find_movie_by_quote():
    """
    Identifies the Oscar-nominated film based on a specific scene description.
    """
    film_title = "The Fugitive"
    year = 1993
    character_name = "Fredrick Sykes (the one-armed man)"
    context = (
        f"In the {year} Oscar-nominated film '{film_title}', the character {character_name} "
        "is trying to escape notice after the murder of Dr. Richard Kimble's wife. "
        "While boarding a crowded bus, he politely says 'Thank you' to the driver. "
        "Dr. Kimble, who is also on the bus, hears his voice and sees him, which "
        "ultimately leads to the film's final confrontation."
    )

    print("The film you're looking for is:")
    print(f"- Film: {film_title}")
    print(f"- Character: {character_name}")
    print("\nHere's the context of the scene:")
    print(context)

# Execute the function to provide the answer
find_movie_by_quote()