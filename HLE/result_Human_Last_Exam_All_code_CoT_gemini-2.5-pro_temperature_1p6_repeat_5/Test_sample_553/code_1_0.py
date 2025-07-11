def find_beatles_song():
    """
    This function identifies the Beatles song based on a specific chord progression.

    The user is asking for a song that moves from the tonic chord to the minor fifth chord.
    - Tonic chord (I): The "home" chord of the key.
    - Minor fifth chord (v): A minor chord built on the fifth degree of the scale.

    An analysis of The Beatles' discography points to one famous example:
    Song: "Michelle"
    Key: F
    Tonic (I) chord: F Major
    Minor fifth (v) chord: C minor

    The progression occurs in the verse:
    - "Michelle, ma belle..." is sung over the F Major (I) chord.
    - "sont des mots qui vont très bien ensemble..." starts on the C minor (v) chord.
    """
    song_title = "Michelle"
    tonic_chord = "F Major"
    minor_fifth_chord = "C minor"

    print(f"The Beatles song that features a prominent chord jump from the tonic to the minor fifth is '{song_title}'.")
    print(f"In this song, the progression is from the tonic chord {tonic_chord} (I) to the minor fifth chord {minor_fifth_chord} (v).")

find_beatles_song()