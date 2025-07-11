def find_beatles_song():
    """
    This function identifies a Beatles song based on a specific chord progression
    and prints the answer.
    """
    # The band in question
    band = "The Beatles"

    # The musical query
    # The progression is from the tonic chord (I) to the minor fifth chord (v).
    # For example, in the key of A major, this would be a move from A major to E minor.
    tonic_chord_roman_numeral = "I"
    minor_fifth_chord_roman_numeral = "v"

    # The answer based on music analysis
    song_title = "Help!"

    # Print the final answer
    print(f"The song by {band} that starts with a chord jump from the tonic to the minor fifth chord is '{song_title}'.")
    print("The progression can be represented in Roman numeral analysis as:")
    print(f"{tonic_chord_roman_numeral} -> {minor_fifth_chord_roman_numeral}")

find_beatles_song()