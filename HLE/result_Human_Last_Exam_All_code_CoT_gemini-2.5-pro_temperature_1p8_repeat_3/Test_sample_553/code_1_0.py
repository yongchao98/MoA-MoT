def find_beatles_song():
    """
    This function identifies a Beatles song based on its opening chord progression.

    The request asks for a song that starts with a jump from the tonic chord (I)
    to the minor fifth chord (v).

    Analysis:
    1. Tonic Chord (I): The 'home' chord of the song's key.
    2. Minor Fifth Chord (v): In a major key, the chord built on the fifth degree
       of the scale is typically major (V). A minor fifth (v) is a less common,
       more harmonically interesting choice.
    3. The Song: "Michelle" from the album Rubber Soul (1965).
    4. The Progression: The intro to "Michelle" is in the key of F.
       - The first chord is F major, which is the tonic (I) of the key of F.
       - The second chord is C minor. C is the fifth note of the F major scale.
         Therefore, C minor is the minor fifth (v) chord.

    The progression is F (I) -> Cmin (v), which matches the user's description perfectly.
    """
    song_title = "Michelle"
    key = "F Major"
    tonic_chord = "F Major (I)"
    minor_fifth_chord = "C minor (v)"

    print(f"The Beatles song that starts with a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print(f"Key of the introduction: {key}")
    print(f"Tonic Chord (I): {tonic_chord.split(' ')[0]}")
    print(f"Minor Fifth Chord (v): {minor_fifth_chord.split(' ')[0]} {minor_fifth_chord.split(' ')[1]}")


find_beatles_song()