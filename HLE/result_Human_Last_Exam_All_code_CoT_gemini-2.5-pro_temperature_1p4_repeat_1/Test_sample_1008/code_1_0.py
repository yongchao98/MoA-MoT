def get_sound_hierarchy():
    """
    This function defines and prints the hierarchy of sound importance in a film's dance scene.
    """
    # The elements are ordered from most to least important.
    sound_elements_in_order = ["Music", "Sound Effects", "Speech"]
    
    # Joining the list into a single string, separated by a comma and space.
    ranked_order = ", ".join(sound_elements_in_order)
    
    print(ranked_order)

get_sound_hierarchy()