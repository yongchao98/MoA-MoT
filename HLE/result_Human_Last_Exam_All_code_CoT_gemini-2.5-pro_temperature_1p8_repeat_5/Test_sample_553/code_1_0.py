def find_beatles_song():
    """
    This function identifies the Beatles song that starts with a jump
    from the tonic chord to the minor fifth chord.
    """
    song_title = "From Me to You"
    artist = "The Beatles"
    key = "C Major"
    
    # The chord progression analysis
    tonic_chord_numeral = "I"
    tonic_chord_name = "C Major"
    
    minor_fifth_chord_numeral = "v7" # The 'v' is the 'number' for the fifth degree, lowercase for minor. The 7th is often included.
    minor_fifth_chord_name = "G minor 7th"

    print(f"The song you're looking for is '{song_title}' by {artist}.")
    print("-" * 30)
    print(f"The song is in the key of {key}.")
    print("The opening chord progression is a jump from the tonic to the minor fifth.")
    print("\nHere is the breakdown of that progression:")
    print(f"Tonic Chord ({tonic_chord_numeral}): {tonic_chord_name}")
    print(f"Minor Fifth Chord ({minor_fifth_chord_numeral}): {minor_fifth_chord_name}")

    print("\nFinal Equation (Harmonic Analysis):")
    print(f"{tonic_chord_numeral} -> {minor_fifth_chord_numeral}")
    print(f"{tonic_chord_name} -> {minor_fifth_chord_name}")

find_beatles_song()
