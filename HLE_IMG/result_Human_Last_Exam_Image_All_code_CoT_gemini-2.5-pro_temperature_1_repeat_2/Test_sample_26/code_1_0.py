def analyze_beethoven_chord():
    """
    Analyzes the harmony of measure 30 from Beethoven's Pathetique Sonata, 1st mvt,
    and provides the correct Roman numeral.
    """
    key = "C minor"
    measure = 30
    
    # Step 1: Identify the notes in the measure
    bass_notes = "C#"
    treble_notes = "E, G#, B"
    all_notes = "C#, E, G#, B"
    
    # Step 2: Determine the chord's identity
    chord_identity = "C-sharp minor seventh (C#m7)"
    
    # Step 3: Analyze the harmonic context and function
    # The chord in m. 30 leads to the D diminished/minor chord in m. 32.
    target_chord_degree = "ii°"
    target_chord_root = "D"
    
    # The chord in m. 30 is built on C#, the leading tone to D.
    function = f"secondary leading-tone chord to the {target_chord_degree}"
    
    # Step 4: Formulate the Roman Numeral
    leading_tone_numeral = "vii°⁷"
    final_roman_numeral = f"{leading_tone_numeral}/{target_chord_degree}"
    
    # Print the analysis
    print("Step-by-step analysis for measure 30 of Beethoven's Pathetique Sonata:")
    print(f"1. The home key of the movement is {key}.")
    print(f"2. The notes present in measure {measure} are {all_notes}.")
    print(f"3. These notes form a {chord_identity} chord.")
    print(f"4. The chord's function is determined by its resolution. It leads to the D minor chord (the '{target_chord_degree}') in measure 32.")
    print(f"5. A chord built on C# that resolves to {target_chord_root} acts as a leading-tone chord.")
    print("6. Therefore, the chord functions as a secondary (or applied) leading-tone chord.")
    print("\nThe correct Roman numeral, which represents the chord's function, is:")
    print(f"vii°⁷ / ii°")

analyze_beethoven_chord()