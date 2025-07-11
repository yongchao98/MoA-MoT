def find_beatles_song():
    """
    Identifies and explains the Beatles song matching the musical criteria.
    """
    song_title = "I'll Follow the Sun"
    band_name = "The Beatles"
    song_key = "C Major"
    tonic_chord = "C Major"
    minor_fifth_chord = "G minor"

    print(f"The song by {band_name} that starts with a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print(f"The song is in the key of {song_key}.")
    print(f"The vocal phrase begins on the tonic chord ({tonic_chord}, the 'I' chord).")
    print(f"It then 'jumps' to the minor fifth chord ({minor_fifth_chord}, the 'v' chord) on the line '...to see I've gone'.")

# Execute the function to print the answer
find_beatles_song()