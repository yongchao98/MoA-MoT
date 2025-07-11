def get_sound_hierarchy():
    """
    Determines and prints the hierarchy of sound importance in a film's dance scene.
    
    In a typical dance scene:
    1. Music is the most important element, as it provides the foundation for choreography and mood.
    2. Sound Effects (like footsteps, cloth sounds) are next, adding realism and texture to the performance.
    3. Speech/Dialogue is usually the least important, often mixed low or omitted entirely to avoid distracting
       from the music and dance.
    """
    
    # The hierarchy of sound importance
    hierarchy = ["Music", "Sound Effects", "Speech"]
    
    # Joining the list into a string separated by a comma and a space
    result = ", ".join(hierarchy)
    
    print(result)

if __name__ == "__main__":
    get_sound_hierarchy()