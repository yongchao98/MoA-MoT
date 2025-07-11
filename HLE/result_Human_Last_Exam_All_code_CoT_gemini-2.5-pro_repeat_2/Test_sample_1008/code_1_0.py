def get_sound_hierarchy():
    """
    Determines and prints the hierarchy of sound importance in film dance scenes.
    """
    # The hierarchy is ordered from most important to least important.
    # 1. Music: The foundation of any dance scene, dictating rhythm and emotion.
    # 2. Sound Effects: Add realism and texture (e.g., footsteps, clothing swishes).
    # 3. Speech: Often minimal or absent to keep the focus on the dance.
    hierarchy = ["Music", "Sound Effects", "Speech"]

    # Format the list for clear, comma-separated output.
    output = ", ".join(hierarchy)
    
    print(output)

if __name__ == "__main__":
    get_sound_hierarchy()