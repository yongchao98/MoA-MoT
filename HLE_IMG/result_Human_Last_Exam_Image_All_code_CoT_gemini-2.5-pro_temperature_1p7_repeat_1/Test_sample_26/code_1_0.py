def analyze_measure_30():
    """
    Analyzes the harmony of measure 30 of Beethoven's "Pathetique" Sonata, 1st mov.
    and determines the correct Roman numeral.
    """
    key = "C minor"
    measure = 30
    
    # Notes identified from the score in measure 30
    notes_in_chord = ["F#", "A", "C", "Eb"]
    
    # Chord analysis
    root = "F#"
    quality = "diminished seventh"
    chord_name = f"{root}°7"
    
    # Functional analysis in the key of C minor
    dominant_of_key = "V (G major)"
    function_of_chord = "The leading-tone chord to the dominant (G)"
    roman_numeral = "vii°⁷/V"
    
    # Print the step-by-step analysis and final answer
    print(f"Analysis of Beethoven 'Pathetique' Sonata, 1st mvt., measure {measure}:")
    print("-" * 60)
    print(f"1. The key of the movement is {key}.")
    print(f"2. The notes present in measure {measure} are: {', '.join(notes_in_chord)}.")
    print(f"3. These notes form an {root} {quality} chord ({chord_name}).")
    print(f"4. In {key}, this chord functions as a secondary leading-tone chord.")
    print(f"   It is the leading-tone chord that resolves to the {dominant_of_key} chord.")
    print("-" * 60)
    print("The correct Roman numeral is therefore:")
    print(f"vii°⁷/V (vii diminished 7 of 5)")
    print("-" * 60)
    
analyze_measure_30()