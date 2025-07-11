def find_beatles_song():
    """
    This function identifies the Beatles song with a specific chord progression
    and prints the answer with a music theory explanation.
    """
    song_title = "Michelle"
    key = "F minor"
    tonic_chord = "F minor (i)"
    minor_fifth_chord = "C minor (v)"

    # The opening guitar passage of "Michelle" is an arpeggiated progression
    # moving from the tonic chord to the minor fifth chord.
    
    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth chord is: {song_title}")
    print("\nExplanation:")
    print(f"The song's verse is in the key of {key}.")
    print(f"The progression jumps from the tonic chord ({tonic_chord}) to the minor fifth chord ({minor_fifth_chord}).")

find_beatles_song()