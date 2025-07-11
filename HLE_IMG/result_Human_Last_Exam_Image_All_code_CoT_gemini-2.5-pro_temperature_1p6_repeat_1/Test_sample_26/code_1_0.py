def analyze_measure_30():
    """
    Analyzes measure 30 of Beethoven's Pathetique Sonata, 1st Movement,
    and provides the correct Roman numeral analysis.
    """
    key = "C minor"
    measure = 30
    
    # Notes found in measure 30
    notes_in_chord = ["F#", "A", "C", "Eb"]
    
    # Analysis of the chord
    chord_name = "F-sharp diminished seventh (F#°7)"
    dominant_of_key = "G major (V)"
    function_of_chord = f"Leading-tone seventh chord to the dominant ({dominant_of_key})"
    roman_numeral = "vii°7/V"
    
    # Print the step-by-step analysis
    print(f"Analysis of Beethoven's Pathetique Sonata, 1st Mvt., Measure {measure}:")
    print(f"1. The key is {key}.")
    print(f"2. The notes in measure {measure} are: {', '.join(notes_in_chord)}.")
    print(f"3. These notes form an {chord_name}.")
    print(f"4. This chord functions as the {function_of_chord}.")
    print(f"5. Therefore, the correct Roman numeral is: {roman_numeral}")

analyze_measure_30()