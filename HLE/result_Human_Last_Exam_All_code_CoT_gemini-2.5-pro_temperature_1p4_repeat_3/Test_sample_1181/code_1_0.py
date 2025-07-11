def find_solaris_character():
    """
    Identifies the character from Tarkovsky's 'Solaris' who misses the sound of leaves.
    This script contains data about the characters and searches for the correct one based on the user's query.
    """
    characters = {
        "Kris": "The protagonist psychologist who, while on the space station, recalls his time on Earth and feels ashamed for not having paid attention to simple things like the sound of rustling leaves.",
        "Hari": "A 'guest' created by the Solaris ocean based on Kris's memories of his deceased wife. She has no memory of Earth herself.",
        "Snaut": "A scientist on the station who has become cynical and weary after his long time in space.",
        "Sartorius": "A cold and pragmatic scientist focused entirely on his research, viewing emotions as a hindrance.",
        "Gibarian": "A friend of Kris who died before his arrival, seen only in recordings where he expresses distress."
    }

    # The detail we are looking for.
    search_term = "rustling leaves"

    # Find the character whose description contains the search term.
    found_character = None
    for name, description in characters.items():
        if search_term in description:
            found_character = name
            break

    if found_character:
        print(f"The character who is ashamed to miss the sound of leaves rustling on Earth is: {found_character}")
    else:
        print("Could not determine the character based on the available information.")

find_solaris_character()