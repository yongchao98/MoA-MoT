def analyze_beethoven_measure_30():
    """
    This function analyzes the harmony of measure 30 from Beethoven's
    Pathetique Sonata, 1st movement, and provides the correct Roman numeral.
    """
    
    # 1. Establish the context
    local_key = "E-flat Major"
    home_key = "C minor"
    
    # 2. Identify the notes in the measure
    notes_in_chord = ["F#", "A", "C", "E-flat"]
    
    # 3. Identify the chord quality
    chord_name = "F-sharp fully-diminished seventh (F#°7)"
    
    # 4. Identify the chord of resolution (in measure 31)
    resolution_chord_root = "G"
    resolution_chord_quality = "minor"
    
    # 5. Determine the function of the resolution chord in the local key
    # In E-flat major (Eb, F, G, Ab, Bb, C, D), G is the 3rd degree (mediant).
    # The diatonic chord on the mediant is G minor.
    resolution_roman_numeral = "iii"
    
    # 6. Determine the function of the chord in measure 30
    # The chord F#°7 is the leading-tone seventh chord of G minor.
    # Therefore, it is a secondary leading-tone chord of the mediant.
    final_roman_numeral = "vii°7/iii"
    
    print("Harmonic Analysis of Beethoven Pathetique Sonata, Mvt. 1, Measure 30:")
    print(f"The local key of this section is {local_key}, the relative major of {home_key}.")
    print(f"The notes in measure 30 are {', '.join(notes_in_chord)}, forming an {chord_name}.")
    print(f"This chord resolves to {resolution_chord_root} {resolution_chord_quality} in the next measure.")
    print(f"In {local_key}, the {resolution_chord_root} {resolution_chord_quality} chord is the mediant, or '{resolution_roman_numeral}'.")
    print(f"The {chord_name} functions as the leading-tone chord to {resolution_chord_root} {resolution_chord_quality}.")
    print("\nThe correct Roman numeral is therefore a secondary leading-tone chord of the mediant.")
    print("\nFinal Answer (numeral breakdown):")
    
    numeral_root = "vii"
    quality_symbol = "°"
    seventh = "7"
    secondary_function_of = "/"
    secondary_target = "iii"
    
    print(f"Leading-tone triad of the secondary key: {numeral_root}")
    print(f"Quality (diminished): {quality_symbol}")
    print(f"Seventh (diminished seventh): {seventh}")
    print(f"Secondary function ('of'): {secondary_function_of}")
    print(f"Target chord (mediant of main key): {secondary_target}")
    
    print(f"\nPutting it all together, the Roman numeral is: {final_roman_numeral}")

analyze_beethoven_measure_30()
