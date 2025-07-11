def find_beatles_song():
    """
    This script identifies and explains the Beatles song that matches
    the user's music theory question.
    """
    song_title = "Things We Said Today"
    artist = "The Beatles"
    key = "A minor"
    tonic_chord = "A minor (Am)"
    minor_fifth_chord = "E minor (Em)"

    # Print the main answer
    print(f"The song by {artist} that starts with a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print("\nHere is the musical explanation:")

    # Print the components of the musical "equation"
    print(f"1. The song's key is {key}.")
    print(f"2. The tonic chord (the 'I' chord) is {tonic_chord}.")
    print(f"3. The minor fifth chord (the 'V' chord in a natural minor scale) is {minor_fifth_chord}.")
    
    print("\nThe song's verse starts by moving directly from Am to Em, fitting the description perfectly.")

if __name__ == "__main__":
    find_beatles_song()