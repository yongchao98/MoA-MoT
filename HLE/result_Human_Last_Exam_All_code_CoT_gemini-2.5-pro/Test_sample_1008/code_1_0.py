def get_sound_hierarchy():
    """
    Determines and prints the hierarchy of sound importance in film dance scenes.
    """
    # The hierarchy is defined with ranks: 1 is most important.
    sound_importance = {
        "Music": 1,
        "Sound Effects": 2,
        "Speech": 3
    }

    # Sort the sound elements based on their importance rank.
    sorted_elements = sorted(sound_importance, key=sound_importance.get)

    # Format the list into the required string format.
    result = ", ".join(sorted_elements)
    
    # Print the final ordered list.
    print(result)

get_sound_hierarchy()