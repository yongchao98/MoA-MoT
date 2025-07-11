def analyze_chord():
    """
    Analyzes the harmony of measure 30 from Beethoven's Pathetique Sonata, 1st mvt.
    and provides the correct Roman numeral.
    """
    key = "G minor"
    notes_in_chord = ["F#", "A", "C", "Eb"]
    chord_root = "F#"
    chord_quality = "fully diminished seventh"
    scale_degree = 7
    roman_numeral_base = "vii"
    roman_numeral_quality = "°7"

    print("Step 1: Identify the local key.")
    print(f"The music has modulated to {key}.\n")

    print("Step 2: Identify the notes in the measure.")
    print(f"The notes in measure 30 are {', '.join(notes_in_chord)}.\n")

    print("Step 3: Determine the chord's root and quality.")
    print(f"These notes form an {chord_root} {chord_quality} chord ({chord_root}{roman_numeral_quality}).\n")

    print("Step 4: Determine the Roman Numeral.")
    print(f"The root of the chord, {chord_root}, is the leading tone (the {scale_degree}th degree) of {key}.")
    print(f"A {chord_quality} chord built on the leading tone is represented by the Roman numeral {roman_numeral_base}{roman_numeral_quality}.\n")

    print("Final Answer:")
    final_numeral = f"{roman_numeral_base}{roman_numeral_quality}"
    print(f"The correct Roman numeral for measure 30 is {final_numeral}.")

analyze_chord()