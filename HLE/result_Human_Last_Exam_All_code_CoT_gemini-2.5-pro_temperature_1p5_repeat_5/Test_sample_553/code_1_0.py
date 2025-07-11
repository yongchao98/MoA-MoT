def find_beatles_song():
    """
    This function identifies the Beatles song based on a specific opening chord progression.
    The progression is a jump from the tonic (I) chord to the minor fifth (v) chord.
    """
    song_title = "You're Going to Lose That Girl"
    key = "G Major"
    tonic_chord = "G Major (I)"
    minor_fifth_chord = "D minor (v)"
    
    print(f"The Beatles song that starts with a chord jump from the tonic to the minor fifth is: {song_title}")
    print(f"The song is in the key of {key}.")
    print(f"The progression is from the tonic chord, {tonic_chord}, to the minor fifth chord, {minor_fifth_chord}.")

find_beatles_song()