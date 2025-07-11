def find_beatles_song_with_progression():
    """
    Identifies and explains the Beatles song that fits the user's
    musical description.
    """
    song_title = "Michelle"
    key = "F minor"
    tonic_chord_name = "F minor"
    minor_fifth_chord_name = "C minor"

    print(f"The song by The Beatles that starts with a jump from the tonic to the minor fifth chord is '{song_title}'.")
    print("\n--- Musical Analysis ---")
    print(f"The song is in the key of {key}.")

    # "Equation" Part 1: The Tonic Chord
    # Using the roman numeral 'i' for the first chord in the progression.
    print(f"1. The tonic (home) chord is {tonic_chord_name}, represented by the Roman numeral 'i'.")

    # "Equation" Part 2: The Minor Fifth Chord
    # Using the roman numeral 'v' for the fifth chord in the progression.
    print(f"2. The minor fifth chord is {minor_fifth_chord_name}, represented by the Roman numeral 'v'.")

    print("\nThis chord progression (i -> v) happens as the main vocal melody begins.")

if __name__ == '__main__':
    find_beatles_song_with_progression()