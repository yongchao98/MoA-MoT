def find_beatles_song():
    """
    This script identifies the Beatles song based on a music theory question.

    The question asks for a song starting with a jump from the tonic to the minor fifth chord.
    - Tonic chord is the first chord of the key (Roman numeral: I).
    - Minor fifth chord is the fifth chord of the key, but played as a minor (Roman numeral: vm).

    The song "Michelle" opens with this exact progression.
    The first line "Michelle, ma belle" is sung over these two chords:
    1. F Major (The tonic, I)
    2. C minor (The minor fifth, vm)
    """

    song_title = "Michelle"
    tonic_chord = "F Major"
    minor_fifth_chord = "C minor"

    print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth chord is '{song_title}'.")
    print(f"The specific progression in the opening is from {tonic_chord} (the tonic, or chord 1) to {minor_fifth_chord} (the minor fifth, or chord 5m).")

# Execute the function to provide the answer
find_beatles_song()