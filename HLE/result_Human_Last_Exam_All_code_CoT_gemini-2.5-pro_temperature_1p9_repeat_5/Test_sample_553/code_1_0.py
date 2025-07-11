def find_beatles_song():
    """
    This function identifies the Beatles song based on a music theory clue
    and prints the analysis.
    """
    song_title = "Girl"
    album = "Rubber Soul"
    key = "C minor"
    tonic_chord = "C minor (the 'i' chord)"
    minor_fifth_chord = "G minor (the 'v' chord)"

    print(f"The Beatles song that starts with a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print(f"It appears on the album '{album}'.")
    print("-" * 20)
    print("Musical Analysis:")
    print(f"Key: {key}")
    print(f"Tonic Chord: {tonic_chord}")
    print(f"Minor Fifth Chord: {minor_fifth_chord}")
    print("The song's distinctive opening guitar part establishes this i-v progression.")

find_beatles_song()