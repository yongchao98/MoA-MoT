def roman_numeral_analysis():
    """
    Provides a step-by-step harmonic analysis of the circled chord
    in measure 8 of Mozart's Fantasy in D minor.
    """
    print("Step 1: Identifying the Key of the Piece")
    key = "D minor"
    print(f"The piece is in {key}.")
    print("-" * 30)

    print("Step 2: Identifying the Notes in the Chord")
    notes = ["E", "G#", "B"]
    print(f"The notes in the circled measure combine to form the pitches: {', '.join(notes)}.")
    print("-" * 30)

    print("Step 3: Determining the Chord Quality")
    chord_name = "E major"
    print(f"The notes {', '.join(notes)} create an {chord_name} triad in root position.")
    print("-" * 30)

    print("Step 4: Analyzing the Chord's Function")
    home_key_v_chord = "A major"
    function_description = "The E major chord strongly leads to the A major chord (the V chord) in the next measure."
    analysis_term = "secondary dominant"
    print(f"In the key of {key}, the dominant (V) chord is {home_key_v_chord}.")
    print(function_description)
    print(f"This makes the E major chord a {analysis_term}.")
    print("-" * 30)

    print("Step 5: Final Roman Numeral")
    roman_numeral = "V/V"
    explanation = f"The E major chord is the dominant of the dominant (V of V). Therefore, the Roman numeral is {roman_numeral}."
    print(explanation)
    print("\nFinal Answer:")
    print(roman_numeral)

roman_numeral_analysis()