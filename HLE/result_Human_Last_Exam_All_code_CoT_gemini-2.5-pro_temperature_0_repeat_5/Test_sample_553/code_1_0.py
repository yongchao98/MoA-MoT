def find_beatles_song():
    """
    Identifies and explains the Beatles song that starts with a
    tonic to minor fifth chord progression.
    """
    song_title = "Michelle"
    artist = "The Beatles"
    key = "F Major"

    # Define the chords based on the key of F Major
    # The Roman numeral for the tonic chord is I
    tonic_chord_numeral = "I"
    tonic_chord_name = "F Major"

    # The Roman numeral for the minor fifth chord is v
    minor_fifth_chord_numeral = "v"
    minor_fifth_chord_name = "C minor"

    print(f"The song by {artist} that features this progression is '{song_title}'.")
    print("\nHere is the musical analysis:")
    print(f"The song is in the key of {key}.")
    print(f"The tonic chord (numeral {tonic_chord_numeral}) is {tonic_chord_name}.")
    print(f"The minor fifth chord (numeral {minor_fifth_chord_numeral}) is {minor_fifth_chord_name}.")
    print("\nIn the opening line, 'Michelle, ma belle', the harmony jumps from the tonic F Major chord to the C minor chord, creating its distinctive, melancholic sound.")

find_beatles_song()