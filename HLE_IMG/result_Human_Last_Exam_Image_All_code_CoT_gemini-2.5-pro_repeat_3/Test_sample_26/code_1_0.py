def analyze_chord():
    """
    Analyzes the chord in measure 30 of Beethoven's Pathetique Sonata, 1st mov.
    and prints the corresponding Roman numeral.
    """
    key = "C minor"
    chord_notes = ["G", "B-natural", "D", "F"]
    bass_note = "G"
    root_note = "G"
    
    # In C minor, G is the 5th scale degree.
    scale_degree = 5
    
    # The chord is a dominant 7th.
    chord_type = 7

    # The triad G-B-D is major, so the Roman numeral is uppercase.
    roman_numeral_base = "V"

    print(f"Analysis of Beethoven's Pathetique Sonata, Op. 13, 1st Mvt., Measure 30:")
    print(f"Key: {key}")
    print(f"Notes in the chord: {', '.join(chord_notes)}")
    print(f"The root of the chord is {root_note}.")
    print(f"The root '{root_note}' is the {scale_degree}th degree of {key}.")
    print(f"The chord quality is a dominant {chord_type}th.")
    print("The bass note is the root, so the chord is in root position.")
    print("-" * 20)
    print("The final Roman numeral is constructed as follows:")
    print(f"Base Roman numeral for scale degree {scale_degree}: {roman_numeral_base}")
    print(f"Add the number for the seventh: {chord_type}")
    print(f"Final Result: {roman_numeral_base}{chord_type}")

analyze_chord()