def find_beatles_song():
    """
    This function identifies and explains the Beatles song
    that matches the user's music theory query.
    """
    song_title = "Got to Get You into My Life"
    key = "G Major"
    tonic_chord = "G Major"
    minor_fifth_chord = "D minor"

    print(f"The Beatles song that starts its verse with a chord jump from the tonic to the minor fifth is: '{song_title}'.")
    print("\n--- Musical Analysis ---")
    print(f"Key: {key}")
    print(f"Tonic Chord (I): {tonic_chord}")
    print(f"Minor Fifth Chord (v): {minor_fifth_chord}")
    
    # Representing the progression as a simple "equation" as requested.
    print("\nChord Equation:")
    print(f"Chord 1 (Tonic) -> Chord 5 (Minor)")
    print(f"The progression is {tonic_chord} -> {minor_fifth_chord}")

find_beatles_song()