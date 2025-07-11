def analyze_mozart_chord():
    """
    Analyzes the circled chord in measure 8 of Mozart's Fantasy in D minor.
    """
    key = "D minor"
    notes_in_chord = ["F#", "A", "C#", "E"]
    chord_name = "F-sharp minor 7th"
    
    following_chord_root = "G"
    following_chord_quality = "minor"
    function_in_key = "subdominant (iv)"
    
    analysis_steps = [
        f"1. The key of the piece is {key}.",
        f"2. The notes in the circled chord are {', '.join(notes_in_chord)}.",
        f"3. Stacking these notes in thirds reveals an {chord_name} chord.",
        f"4. This chord resolves to a {following_chord_root} {following_chord_quality} chord in the next measure, which is the {function_in_key} of {key}.",
        f"5. The root of the chord, F#, is the leading tone to the next chord's root, {following_chord_root}. This makes its function a 'secondary leading-tone chord'.",
        "6. We construct the Roman Numeral as follows:",
        "   - The root is the leading tone of the chord of resolution -> 'vii'",
        "   - The chord quality is minor with a seventh -> lowercase 'vii' with a '7' -> 'vii⁷'",
        "   - The chord resolves to the subdominant (iv) -> '/iv'",
        f"7. Combining these gives the final Roman numeral."
    ]
    
    print("Step-by-step analysis of the chord:")
    for step in analysis_steps:
        print(step)
        
    final_numeral = "vii⁷/iv"
    print("\nFinal Roman Numeral Equation:")
    print(f"({chord_name}) in the key of ({key}) functions as a secondary leading-tone chord to ({function_in_key}), so the Roman Numeral is {final_numeral}")

analyze_mozart_chord()