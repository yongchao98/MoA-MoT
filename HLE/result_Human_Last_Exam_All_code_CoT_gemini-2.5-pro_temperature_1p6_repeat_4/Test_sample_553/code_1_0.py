import sys

def find_beatles_song():
    """
    This function analyzes a music theory question about The Beatles' compositions
    and provides the answer.

    The question: Which song by The Beatles starts with a chord jump from the
    tonic to the minor fifth chord?

    Analysis:
    - Tonic chord is denoted by the Roman numeral: I
    - Minor fifth chord is denoted by the Roman numeral: v
    - The progression is I -> v. This is an unusual progression that often involves
      borrowing a chord from the parallel minor key.
    - The song "Michelle" contains the most famous example of this. The bridge,
      in F Major (the tonic, I), transitions to a C minor chord (the minor fifth, v).

    The equation representing the chord jump is: Tonic (I) -> minor fifth (v)
    """

    song_name = "Michelle"
    tonic_chord_numeral = "I"
    minor_fifth_chord_numeral = "v"

    # Print the answer to the console.
    print(f"The Beatles song featuring a prominent chord jump from the tonic ({tonic_chord_numeral}) to the minor fifth ({minor_fifth_chord_numeral}) is '{song_name}'.")
    # The "equation" requested refers to the chord progression numerals.
    print(f"The chord progression 'equation' is: {tonic_chord_numeral} -> {minor_fifth_chord_numeral}")

# Execute the function to find and print the answer.
find_beatles_song()
