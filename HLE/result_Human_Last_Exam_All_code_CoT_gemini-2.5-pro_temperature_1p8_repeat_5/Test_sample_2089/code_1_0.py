def find_movie_scene():
    """
    This function provides the answer to the user's movie trivia question.
    """
    movie_title = "Carlito's Way"
    scene_description = (
        "In the Oscar-nominated film '{title}', the character Lalin reveals he is a "
        "government informant by saying 'Thank you' into his hidden microphone as he "
        "boards a bus. This act seals his fate in the eyes of the protagonist, Carlito Brigante."
    ).format(title=movie_title)
    
    print(scene_description)

find_movie_scene()