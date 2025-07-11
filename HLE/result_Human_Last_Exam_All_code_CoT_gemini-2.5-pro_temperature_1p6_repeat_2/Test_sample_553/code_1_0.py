def find_beatles_song():
    """
    This function analyzes the user's music theory question and provides the answer.
    """
    song_title = "I'll Be Back"
    album = "A Hard Day's Night"
    key = "A Major"
    intro_key = "A Minor"
    
    # Chord Analysis
    tonic_chord_of_intro = "A Minor (Am)"
    minor_fifth_chord = "E Minor (Em)"
    
    # Explanation
    print("Analyzing the chord progression: Tonic (I) to minor fifth (v).")
    print("This is an uncommon progression, often involving a shift between major and minor modes.")
    print("-" * 20)
    print(f"The song that fits this description is '{song_title}' from the album '{album}'.")
    print("-" * 20)
    print("Here is the step-by-step analysis:")
    print(f"1. The overall key of the song is {key}.")
    print(f"2. However, the introduction is played in the parallel minor key, which is {intro_key}.")
    print(f"3. In {intro_key}, the tonic chord (i) is {tonic_chord_of_intro}.")
    print(f"4. The fifth note of the A minor scale is E. The chord built on this note is {minor_fifth_chord}.")
    print(f"5. This chord, {minor_fifth_chord}, is the minor fifth (v) relative to the {intro_key} tonic.")
    print("-" * 20)
    print("The intro progression starts with the following chord change:")
    print(f"{tonic_chord_of_intro} -> {minor_fifth_chord}")
    print("This is a clear example of a song starting with a chord jump from the tonic to the minor fifth.")
    print("-" * 20)
    print("Final Answer:")
    print(song_title)

find_beatles_song()