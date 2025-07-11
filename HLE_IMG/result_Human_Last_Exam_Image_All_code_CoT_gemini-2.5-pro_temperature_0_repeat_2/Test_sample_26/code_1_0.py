def analyze_chord():
    """
    Analyzes the chord in measure 30 of Beethoven's Pathetique Sonata, 1st mvt.
    and determines its Roman numeral in the key of C minor.
    """
    key = "C minor"
    measure_number = 30
    
    # Notes identified from the score in measure 30
    notes_in_chord = ["G", "B-natural", "D", "F-sharp"]
    
    # Step-by-step analysis
    print(f"Analyzing Measure {measure_number} of Beethoven's Pathetique Sonata...")
    print(f"The key of the piece is {key}.")
    print(f"The notes in the chord are: {', '.join(notes_in_chord)}.")
    
    root = "G"
    scale_degree = 5
    print(f"The root of the chord is {root}.")
    print(f"In the key of {key}, {root} is the {scale_degree}th scale degree, which is the Dominant.")
    
    # The triad G-B(natural)-D is a major triad.
    # In Roman numeral analysis, a major chord on the dominant is represented by an uppercase 'V'.
    roman_numeral = "V"
    
    print(f"The core triad is G - B-natural - D, which is a major triad.")
    print("A major chord built on the 5th scale degree is represented by the uppercase Roman numeral 'V'.")
    print("The F-sharp is a chromatic alteration (a major seventh) but does not change the chord's fundamental function as the dominant.")
    print("\nTherefore, the correct Roman numeral is:")
    print(roman_numeral)

analyze_chord()