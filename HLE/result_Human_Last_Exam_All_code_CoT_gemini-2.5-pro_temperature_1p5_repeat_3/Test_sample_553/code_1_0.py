def find_beatles_song():
    """
    Identifies the Beatles song that starts with a jump from the tonic
    to the minor fifth chord.
    """
    song_title = "Michelle"
    key = "F Major"
    tonic_chord = "F Major (I)"
    minor_fifth_chord = "C minor (v)"

    print(f"The Beatles song that starts with a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print(f"In the key of {key}, the progression is a jump from the tonic {tonic_chord} to the minor fifth {minor_fifth_chord}.")

find_beatles_song()