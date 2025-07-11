def find_beatles_song():
    """
    This script identifies and explains the Beatles song
    that fits the user's description.
    """
    song_name = "Michelle"
    album_name = "Rubber Soul"
    key = "F"
    tonic_chord = "F Major (I)"
    minor_fifth_chord = "C minor (v)"

    explanation = f"""
The song by The Beatles that starts with a chord jump from the tonic to the minor fifth chord is '{song_name}'.

Here is a step-by-step explanation:

1.  The song is '{song_name}' from the album '{album_name}'.
2.  The beautiful instrumental introduction establishes the musical key of '{key}' as the tonic center (the 'home' chord).
3.  When the vocals begin ("Michelle, ma belle..."), the first chord of the verse is a {minor_fifth_chord}.
4.  This creates a beautiful and musically interesting jump from the established {tonic_chord} to the {minor_fifth_chord}.
"""
    print(explanation)

find_beatles_song()