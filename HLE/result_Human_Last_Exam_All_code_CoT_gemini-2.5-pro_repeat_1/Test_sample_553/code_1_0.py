def solve_beatles_riddle():
    """
    Identifies and explains the Beatles song that uses a
    tonic to minor fifth chord progression.
    """
    song_title = "Michelle"
    key_of_section = "F minor"
    tonic_chord = "F minor"
    roman_numeral_tonic = "i"
    minor_fifth_chord = "C minor"
    roman_numeral_minor_fifth = "v"

    print(f"The song by The Beatles that features a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print("\nThis progression is a key feature of the song's bridge ('I love you, I love you, I love you...').")
    print(f"In the key of {key_of_section} for that section, the progression is as follows:")
    print("-" * 30)
    print(f"The Tonic Chord ({roman_numeral_tonic}): {tonic_chord}")
    print(f"The Minor Fifth Chord ({roman_numeral_minor_fifth}): {minor_fifth_chord}")
    print("-" * 30)
    print(f"The chord change from '{tonic_chord}' to '{minor_fifth_chord}' on the words 'I love' is the moment you are asking about.")

solve_beatles_riddle()