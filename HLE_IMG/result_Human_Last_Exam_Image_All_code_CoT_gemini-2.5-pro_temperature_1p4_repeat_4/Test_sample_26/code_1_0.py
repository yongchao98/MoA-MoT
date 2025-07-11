def analyze_chord():
    """
    Analyzes the chord in measure 30 of Beethoven's Pathetique Sonata, Mvt. 1.
    """
    # Step 1: Define the musical context
    key = "E minor"
    measure = 30
    print(f"Step 1: The analysis is performed in the context of the temporary key, which is {key}.")

    # Step 2: Identify the notes in the chord
    # Treble clef has C# and E. Bass clef has F# and A.
    notes = ["F#", "A", "C#", "E"]
    print(f"Step 2: The notes in measure {measure} are {', '.join(notes)}.")

    # Step 3: Determine the chord's root and quality
    # Stacking the notes in thirds from F#: F# - A - C# - E
    root = "F#"
    # F# to A is a minor third. F# to C# is a perfect fifth. -> F# minor triad.
    triad_quality = "minor"
    # F# to E is a minor seventh.
    seventh_quality = "minor seventh"
    chord_name = f"{root} {triad_quality} seventh ({root}m7)"
    print(f"Step 3: The chord is an {chord_name}. The root is {root}.")

    # Step 4: Determine the harmonic function
    # In the key of E minor, F# is the second scale degree (supertonic).
    scale_degree = "ii"
    print(f"Step 4: The root ({root}) is the supertonic (second scale degree) of {key}.")

    # Step 5: Construct the Roman numeral
    # - Lowercase 'ii' because the triad is minor.
    # - '7' is added for the minor seventh.
    # - The chord is in root position (F# is in the bass), so no inversion figures are needed.
    # This chord (iim7) is a common pre-dominant chord, derived from the melodic minor scale, leading to the V7 (B7) in the next measure.
    roman_numeral = "ii⁷"
    print(f"Step 5: The Roman numeral is '{scale_degree}' because the triad is minor. A '7' is added for the seventh.")
    print("-" * 20)
    print(f"The final Roman numeral for the chord in measure {measure} is {roman_numeral}.")

analyze_chord()