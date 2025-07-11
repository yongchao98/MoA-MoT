def find_beatles_song():
    """
    This function identifies the Beatles song that starts with a chord jump
    from the tonic to the minor fifth and explains the music theory behind it.
    """
    song_title = "Michelle"
    key = "F Major"
    tonic_chord = "F Major (I)"
    minor_fifth_chord = "C minor (v)"

    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print("\nHere is the breakdown:")
    print(f"1. The song is in the key of {key}.")
    print(f"2. The first chord, the tonic (I), is {tonic_chord.split(' ')[0]}.")
    print(f"3. The second chord is {minor_fifth_chord.split(' ')[0]}. In the key of {key}, the fifth note is C. A C minor chord is therefore the minor fifth (v).")
    print(f"\nThe progression is {tonic_chord} -> {minor_fifth_chord}.")

find_beatles_song()