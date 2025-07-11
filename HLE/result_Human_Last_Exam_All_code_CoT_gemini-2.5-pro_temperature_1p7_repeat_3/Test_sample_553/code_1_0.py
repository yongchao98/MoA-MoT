def solve_beatles_song_query():
    """
    This script identifies the Beatles song based on a specific
    opening chord progression and explains the music theory behind it.
    """
    song_title = "Things We Said Today"
    key = "A minor"
    
    # Define the components of the chord progression "equation"
    tonic_chord = {
        "roman_numeral": "i",
        "chord_name": "A minor"
    }
    
    minor_fifth_chord = {
        "roman_numeral": "v",
        "chord_name": "E minor"
    }

    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print(f"\nThe song's opening is in the key of {key}.")
    print("\nThe musical 'equation' for this opening progression is:")

    # Printing each component of the 'equation' as requested.
    print(f"Term 1 (Tonic Chord): {tonic_chord['roman_numeral']} which is {tonic_chord['chord_name']}")
    print(f"Term 2 (Minor Fifth Chord): {minor_fifth_chord['roman_numeral']} which is {minor_fifth_chord['chord_name']}")

solve_beatles_song_query()