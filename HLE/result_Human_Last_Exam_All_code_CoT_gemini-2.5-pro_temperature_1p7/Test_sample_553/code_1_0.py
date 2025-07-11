def find_beatles_song():
    """
    Identifies and prints the Beatles song that starts with a
    tonic to minor fifth chord progression.
    """
    song_title = "Magical Mystery Tour"
    key = "E Major"
    tonic_chord = "E Major (the tonic, or I chord)"
    minor_fifth_chord = "B minor (the minor fifth, or vm chord)"

    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print("\nThe song is in the key of E Major.")
    print("The chord 'equation' or progression is:")
    print(f"1. Tonic Chord: {tonic_chord}")
    print(f"2. Minor Fifth Chord: {minor_fifth_chord}")

if __name__ == "__main__":
    find_beatles_song()