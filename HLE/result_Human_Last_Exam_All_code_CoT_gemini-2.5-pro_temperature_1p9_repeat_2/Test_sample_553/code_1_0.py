def find_beatles_song():
    """
    Identifies and explains the Beatles song that starts with a specific chord progression.
    """
    song_title = "Michelle"
    key_context = "F minor"
    tonic_chord = "F minor"
    minor_fifth_chord = "C minor"

    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth chord is '{song_title}'.")
    print("\nHere is a brief musical analysis of the song's introduction:")
    print(f"The song's intro is in the key of {key_context}.")
    print(f"It begins on the tonic chord: {tonic_chord} (represented as 'i').")
    print(f"It then jumps to the minor fifth chord: {minor_fifth_chord} (represented as 'v').")
    print("\nThis creates the distinctive, melancholic opening that the song is known for.")

find_beatles_song()