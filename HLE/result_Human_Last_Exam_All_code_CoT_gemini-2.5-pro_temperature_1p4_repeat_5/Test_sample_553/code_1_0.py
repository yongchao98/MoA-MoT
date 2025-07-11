def find_beatles_song():
    """
    Identifies the Beatles song that starts with a I -> v chord progression
    and prints the details.
    """
    song_title = "Got to Get You Into My Life"
    key = "G Major"
    
    # In the key of G Major:
    # The tonic chord (I) is G Major. The Roman numeral is I, representing the 1st degree.
    tonic_chord_name = "G Major"
    tonic_degree = 1
    
    # The minor fifth chord (v) is a minor chord built on the 5th degree of the scale.
    # The 5th degree is D, so the chord is D minor.
    minor_fifth_chord_name = "D minor"
    minor_fifth_degree = 5

    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth is: '{song_title}'")
    print(f"This progression happens in the opening brass fanfare of the song.")
    print("-" * 40)
    print(f"The song's key: {key}")
    print("The musical 'equation' for the chord progression is:")
    
    # The final output needs to include the numbers from the 'equation'
    print(f"Chord 1 (Tonic): {tonic_chord_name} (Degree: {tonic_degree})")
    print(f"Chord 2 (Minor Fifth): {minor_fifth_chord_name} (Degree: {minor_fifth_degree}, minor quality)")
    
find_beatles_song()