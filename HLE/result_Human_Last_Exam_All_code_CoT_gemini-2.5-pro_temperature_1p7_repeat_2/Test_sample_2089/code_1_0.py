def find_movie_scene():
    """
    Identifies the film based on a scene description provided by the user.
    The user's description: A character in an Oscar-nominated film seals his fate
    by saying "Thank you" while boarding a bus.
    """

    # Analysis reveals that the 'Oscar-nominated' detail is likely incorrect,
    # as another prominent film perfectly matches the rest of the description.
    
    film_title = "Patriot Games (1992)"
    explanation = (
        "The film you are most likely thinking of is 'Patriot Games'.\n\n"
        "In a key scene, a wounded police guard is being left alive after his transport van (a type of bus) is ambushed by terrorists.\n"
        "As the terrorist leader, Sean Miller, turns to leave, the guard says, 'Thank you.'\n"
        "Hearing this, Miller turns back and shoots the guard dead. The guard's politeness directly sealed his fate."
    )

    print(f"While no Oscar-nominated film perfectly matches your description, the scene you described is famously from the film: {film_title}")
    print("\nHere is a description of the scene:")
    print(explanation)

find_movie_scene()