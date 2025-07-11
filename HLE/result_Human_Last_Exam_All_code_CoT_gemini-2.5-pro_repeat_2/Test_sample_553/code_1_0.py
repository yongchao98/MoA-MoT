def find_beatles_song():
    """
    This function identifies the Beatles song based on a specific chord progression.
    The progression is from the tonic chord to the minor fifth chord.

    Analysis:
    - The song "And I Love Her" has a verse section in the key of C# minor.
    - The tonic chord (i) is C# minor.
    - The fifth degree of the C# minor scale is G#. A minor chord built on this is G# minor (v).
    - The song's verse starts by alternating directly between C# minor and G# minor.
    - This perfectly matches the description "a chord jump from the tonic to the minor fifth chord".
    """
    song_title = "And I Love Her"
    tonic_chord = "C# minor"
    minor_fifth_chord = "G# minor"

    print(f"The song by The Beatles that fits this description is: {song_title}")
    print(f"The opening of the verse features a chord jump from the tonic ({tonic_chord}) to the minor fifth ({minor_fifth_chord}).")

find_beatles_song()