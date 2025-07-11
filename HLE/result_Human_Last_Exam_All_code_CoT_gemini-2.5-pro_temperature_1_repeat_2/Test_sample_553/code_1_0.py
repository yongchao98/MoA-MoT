def find_beatles_song():
    """
    This function identifies and explains the Beatles song that matches
    the user's music theory query.
    """
    band = "The Beatles"
    song_title = "And I Love Her"
    description = "a chord jump from the tonic to the minor fifth chord"

    # Musical Analysis details
    key = "C# minor"
    tonic_chord_symbol = "i"
    tonic_chord_name = "C# minor"
    minor_fifth_chord_symbol = "v"
    minor_fifth_chord_name = "G# minor"

    print(f"The song by {band} that starts with {description} is '{song_title}'.")
    print("\n--- Musical Analysis ---")
    print(f"The song's famous opening guitar riff establishes the key of {key}.")
    print("The two chords that make up this iconic intro are:")

    # Using the prompt's instruction to output each part of the "equation" or progression
    print(f"Chord 1 (The Tonic): The '{tonic_chord_symbol}' chord, which is {tonic_chord_name}.")
    print(f"Chord 2 (The Minor Fifth): The '{minor_fifth_chord_symbol}' chord, which is {minor_fifth_chord_name}.")

    print(f"\nThis progression from {tonic_chord_symbol} to {minor_fifth_chord_symbol} ({tonic_chord_name} -> {minor_fifth_chord_name}) perfectly fits the description.")

find_beatles_song()