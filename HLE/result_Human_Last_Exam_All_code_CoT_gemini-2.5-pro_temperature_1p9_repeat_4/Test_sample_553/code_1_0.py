def find_beatles_song():
    """
    This function identifies the Beatles song based on a specific harmonic progression.
    The progression is from the tonic chord to the minor fifth chord (i -> v).

    Analysis:
    - The song is "I'll Be Back".
    - The key is A minor.
    - The tonic chord (i) is A minor.
    - The minor fifth chord (v) is E minor.
    - The song starts with this exact progression.
    """
    song_title = "I'll Be Back"
    key = "A minor"
    tonic_chord = "A minor (i)"
    minor_fifth_chord = "E minor (v)"

    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth chord is '{song_title}'.")
    print(f"The song is in the key of {key}.")
    print(f"It starts with a progression from the tonic chord, {tonic_chord}, to the minor fifth chord, {minor_fifth_chord}.")

find_beatles_song()