def find_film_and_scene():
    """
    This function identifies the film and character based on the user's query.
    """
    film_title = "The Shawshank Redemption"
    character_name = "Brooks Hatlen"
    scene_context = (
        "After being paroled from prison after 50 years, "
        "is overwhelmed by the modern world. While boarding a bus, he politely says 'Thank you' "
        "to the driver. This small act highlights his deep-seated institutionalization and inability to cope "
        "with his newfound freedom, which ultimately leads to his tragic fate."
    )

    print(f"The Oscar-nominated film is '{film_title}'.")
    print(f"The character is {character_name}, who seals his fate by saying \"Thank you\" while boarding a bus.")
    print(f"Context: The character, {character_name}, {scene_context}")

find_film_and_scene()