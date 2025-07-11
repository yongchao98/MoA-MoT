def find_beatles_song():
    """
    This function identifies and explains the Beatles song
    that begins with a tonic to minor fifth chord progression.
    """
    song_title = "Michelle"
    key = "F minor"
    tonic_chord_name = "F minor"
    tonic_chord_roman = "i"
    minor_fifth_chord_name = "C minor"
    minor_fifth_chord_roman = "v"

    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth chord is '{song_title}'.")
    print("-" * 40)
    print("Here is the breakdown of the musical progression:")
    print(f"The key of the song is: {key}")
    print("\nIn the final equation, the opening chords are:")
    # Using '->' to represent the jump from one chord to the next.
    print(f"{tonic_chord_name} ({tonic_chord_roman}) -> {minor_fifth_chord_name} ({minor_fifth_chord_roman})")

find_beatles_song()