def identify_beatles_song():
    """
    Identifies the Beatles song that starts with a chord jump
    from the tonic to the minor fifth chord.
    """
    song_title = "Things We Said Today"
    key = "A minor"
    
    # In the key of A minor, the chords are derived from the A natural minor scale (A-B-C-D-E-F-G).
    # The tonic (first degree) chord is A minor.
    # The fifth degree chord (built on E) is E minor.
    tonic_chord_name = "A minor"
    tonic_chord_symbol = "Am"
    
    minor_fifth_chord_name = "E minor"
    minor_fifth_chord_symbol = "Em"
    
    print(f"The Beatles song that starts with a chord jump from the tonic to the minor fifth is: '{song_title}'")
    print("-" * 50)
    print(f"The song is in the key of {key}.")
    print("The progression in question occurs at the start of the verse.")
    print("\nHere is the final breakdown of the chord 'equation':")
    print(f"1. The Tonic Chord: {tonic_chord_name} ({tonic_chord_symbol})")
    print(f"2. The Minor Fifth Chord: {minor_fifth_chord_name} ({minor_fifth_chord_symbol})")

identify_beatles_song()