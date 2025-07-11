def find_beatles_song():
    """
    This function identifies and explains the Beatles song that fits the user's musical description.
    This is a knowledge-based task, so the information is stored directly in the code.
    """
    song_title = "Michelle"
    explanation = (
        "The Beatles song that features a prominent chord jump from the tonic to a minor fifth chord is 'Michelle'.\n\n"
        "The musical term for this is called modal interchange. The song is in the key of F. "
        "The 'tonic' chord is F major.\n"
        "The standard 'fifth' chord should be C major.\n"
        "However, in the section with the lyrics 'I want you, I want you, I want you...', "
        "the progression moves from the tonic F major to a C minor chord (the minor fifth), "
        "which creates its distinctive and melancholy sound."
    )

    print(f"The song is: {song_title}")
    print("\n--- Explanation ---")
    print(explanation)

find_beatles_song()