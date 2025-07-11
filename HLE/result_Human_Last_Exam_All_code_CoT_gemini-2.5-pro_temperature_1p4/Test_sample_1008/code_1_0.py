def get_sound_hierarchy():
    """
    Defines and prints the hierarchy of sound importance in a film's dance scene.
    """
    # The elements are ordered from most important to least important.
    hierarchy = ["Music", "Sound Effects", "Speech"]
    
    # The final ranked list, formatted as requested.
    print("The hierarchy of sound importance in a dance scene is:")
    print(", ".join(hierarchy))

get_sound_hierarchy()