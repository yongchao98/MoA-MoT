def find_beatles_song():
    """
    This function analyzes and identifies the Beatles song that starts with a
    progression from the tonic to the minor fifth chord.
    """
    song = "Michelle"
    key = "F Major"
    tonic_chord = "F Major (I)"
    minor_fifth_chord = "C minor (v)"

    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth chord is: {song}")
    print("\nHere is the musical analysis:")
    print(f"Key: {key}")
    print(f"The first chord, the tonic (I), is: {tonic_chord.split(' ')[0]}")
    print(f"The second chord, the minor fifth (v), is: {minor_fifth_chord.split(' ')[0]}")
    print(f"\nThe progression is {tonic_chord} -> {minor_fifth_chord}")

find_beatles_song()