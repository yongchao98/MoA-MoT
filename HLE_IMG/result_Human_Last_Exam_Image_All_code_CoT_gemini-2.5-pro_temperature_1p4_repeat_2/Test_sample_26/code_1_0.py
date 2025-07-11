def analyze_chord():
    """
    Analyzes the harmony of measure 30 of Beethoven's "Pathetique" Sonata, mvt. 1.
    """
    # 1. Notes in the measure
    root = "F#"
    third = "A"
    fifth = "C"
    seventh = "Eb"

    # 2. Chord Identification
    chord_name = "F# fully diminished seventh (F#°7)"

    # 3. Harmonic Function Analysis
    # The root F# is the leading tone to G.
    # The key of the piece is C minor. The dominant (V) chord is G minor.
    # Therefore, F#°7 is the leading-tone chord of the dominant.
    analysis = "The chord is a " + chord_name + "." \
        " Its root, " + root + ", is the leading tone to G." \
        " In the home key of C minor, G is the dominant (V)." \
        " This makes the chord the leading-tone seventh of the dominant."

    # 4. Roman Numeral Construction
    leading_tone_numeral = "vii"
    quality_and_inversion = "°7"
    applied_to = "/V"
    final_roman_numeral = leading_tone_numeral + quality_and_inversion + applied_to

    # Print the result
    print(f"The notes in measure 30 are {root}, {third}, {fifth}, and {seventh}.")
    print(f"These notes form an {chord_name}.")
    print(f"Functionally, this chord is the leading-tone chord to G minor.")
    print(f"In the context of the overall key (C minor), G minor is the dominant (V).")
    print(f"Therefore, the correct Roman numeral is {leading_tone_numeral}{quality_and_inversion}{applied_to}.")


if __name__ == "__main__":
    analyze_chord()
