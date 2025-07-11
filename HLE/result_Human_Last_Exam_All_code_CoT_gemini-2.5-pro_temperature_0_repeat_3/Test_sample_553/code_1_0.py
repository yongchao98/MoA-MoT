def find_beatles_song():
    """
    This function identifies and provides details about the Beatles song
    that starts with a chord jump from the tonic to the minor fifth.
    """
    song_title = "Michelle"
    album = "Rubber Soul"
    year = 1965
    key = "F Major"
    
    # In the key of F Major:
    # The tonic chord (I) is F Major.
    # The fifth scale degree is C. A normal V chord would be C Major.
    # A minor fifth chord (v) is C minor.
    
    tonic_chord_name = "F Major"
    tonic_chord_number = 1 # Represents the 'I' chord
    
    minor_fifth_chord_name = "C minor"
    minor_fifth_chord_number = 5 # Represents the 'v' chord
    
    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print("\n--- Musical Details ---")
    print(f"Album: {album} ({year})")
    print(f"Key of the song: {key}")
    print("\nThe opening chord progression is:")
    
    # The prompt requested to output each number in the final equation.
    # We will represent the chord progression as an equation of chord numbers.
    print(f"Chord {tonic_chord_number} ({tonic_chord_name}) -> Chord {minor_fifth_chord_number} ({minor_fifth_chord_name})")
    print("\nThis progression from the tonic (I) to a minor fifth (v) is a key element of the song's unique character.")

find_beatles_song()