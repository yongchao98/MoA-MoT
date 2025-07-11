def find_beatles_song():
    """
    This function identifies the Beatles song based on a specific chord progression.

    The question asks for a song starting with a jump from the tonic (I) chord
    to the minor fifth (v) chord.

    Analysis:
    - Song: "Michelle"
    - Key (Chorus): F Major
    - Tonic Chord (I): F Major
    - Fifth Degree of F: C
    - Minor Fifth Chord (v): C minor

    The chorus of "Michelle" prominently features the progression F Major -> C minor,
    which is a classic example of a I -> v chord change.
    """
    song_title = "Michelle"
    
    # The progression in the key of F Major is F (I) to Cm (v).
    key = "F Major"
    tonic_chord = "F Major (I)"
    minor_fifth_chord = "C minor (v)"

    print(f"The song by The Beatles that features a prominent chord jump from the tonic to the minor fifth is: {song_title}")
    print(f"In the chorus, the progression goes from {tonic_chord} to {minor_fifth_chord}.")

find_beatles_song()