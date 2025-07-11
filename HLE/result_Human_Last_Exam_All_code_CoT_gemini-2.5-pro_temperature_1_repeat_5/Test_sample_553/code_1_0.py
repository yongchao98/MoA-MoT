def find_beatles_song():
    """
    This function identifies the Beatles song based on a specific chord progression
    and prints the analysis.
    """
    
    # Song and key information
    song_title = "I'll Be Back"
    album = "A Hard Day's Night"
    key = "A Major"
    
    # Analysis of the chord progression described
    tonic_chord_name = "A Major"
    tonic_chord_numeral = "I"
    
    minor_fifth_chord_name = "E minor"
    minor_fifth_chord_numeral = "v"
    
    # Print the final answer and explanation
    print(f"The Beatles song you're looking for is: '{song_title}' from the album '{album}'.")
    print("\n--- Analysis of the Chord Progression ---")
    print(f"The song is in the key of {key}.")
    print(f"The 'tonic' chord is {tonic_chord_name} ({tonic_chord_numeral}).")
    print(f"The 'minor fifth chord' is {minor_fifth_chord_name} ({minor_fifth_chord_numeral}).")
    
    # As requested, here is the progression represented as an "equation"
    print("\nThe final progression is:")
    print(f"{tonic_chord_numeral} -> {minor_fifth_chord_numeral}")
    print(f"({tonic_chord_name} -> {minor_fifth_chord_name})")

    print("\nNote: While this progression is a key feature of the song, the intro famously begins with a different major-to-minor shift (A Major to A minor).")

# Execute the function to provide the answer
find_beatles_song()