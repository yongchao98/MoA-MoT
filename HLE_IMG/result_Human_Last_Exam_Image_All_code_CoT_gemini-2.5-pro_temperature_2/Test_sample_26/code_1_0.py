def analyze_chord():
    """
    Analyzes the chord in measure 30 of Beethoven's "Pathetique" Sonata, 1st mov.
    and determines its Roman numeral.
    """
    key = "G minor"
    bass_notes = ["D"]
    treble_notes = ["F#", "A", "C"]
    
    # Combine notes and find the unique set
    all_notes = sorted(list(set(bass_notes + treble_notes)))
    
    root = "D"
    third = "F#"
    fifth = "A"
    seventh = "C"
    
    chord_name = "D dominant 7th (D7)"
    
    # G minor scale: G, A, Bb, C, D, Eb, F
    # Dominant is the 5th degree.
    scale_degree = 5
    
    # In minor keys, the dominant chord is typically major, using the raised 7th
    # of the harmonic or melodic minor scale as its third.
    # The G harmonic minor scale is G, A, Bb, C, D, Eb, F#.
    # The chord D-F#-A is the major dominant (V).
    # The added C is the minor seventh.
    
    roman_numeral_base = "V"
    roman_numeral_quality = "7"
    final_roman_numeral = roman_numeral_base + roman_numeral_quality
    
    print("Step-by-step analysis for measure 30:")
    print(f"1. The key signature indicates the key of {key}.")
    print(f"2. The notes present in the measure are {', '.join(all_notes)}.")
    print(f"3. These notes form a {chord_name} chord.")
    print(f"4. The root of the chord is {root}.")
    print(f"5. In the key of {key}, {root} is the {scale_degree}th scale degree (the dominant).")
    print(f"6. The chord is a major dominant (V) with a minor seventh (7), because the leading tone (F#) is raised.")
    print(f"\nTherefore, the correct Roman numeral is: {final_roman_numeral}")

analyze_chord()