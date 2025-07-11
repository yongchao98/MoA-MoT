def solve_beatles_riddle():
    """
    This function identifies the Beatles song based on a music theory description
    and prints the analysis.
    """
    song_title = "Help!"
    key = "A Major"
    tonic_chord_name = "A Major"
    tonic_chord_numeral = "I"
    following_chord_name = "F# minor"
    following_chord_numeral = "vi"
    actual_minor_fifth_name = "E minor"
    actual_minor_fifth_numeral = "v"

    print("The Beatles song that famously starts with a chord jump from the tonic to a minor harmony is 'Help!'.")
    print("\nHere is the musical breakdown:")
    print(f"The song is in the key of {key}.")
    print(f"It starts with the first word, 'Help!', sung over the powerful tonic chord: {tonic_chord_name} (the '{tonic_chord_numeral}' chord).")
    print(f"The harmony then makes a dramatic 'jump' to a descending progression, with the first main change being to {following_chord_name} (the '{following_chord_numeral}' chord).")
    print(f"This immediate shift from the strong major tonic to a sadder-sounding minor chord is a defining feature of the song's opening.")
    print(f"\nNote: While the jump is not to the technical 'minor fifth' chord ({actual_minor_fifth_name}, the '{actual_minor_fifth_numeral}' chord), this is the most famous example that fits the description.")

solve_beatles_riddle()