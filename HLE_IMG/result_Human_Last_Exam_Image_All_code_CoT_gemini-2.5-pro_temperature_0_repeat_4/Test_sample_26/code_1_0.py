def analyze_measure():
    """
    Analyzes measure 30 of Beethoven's "Pathetique" Sonata, 1st movement,
    to determine the correct Roman numeral.
    """
    key = "C minor"
    measure_number = 30
    
    # Step 1: Identify the key
    print(f"Step 1: The key signature has three flats (Bb, Eb, Ab), indicating the key of {key}.")
    
    # Step 2: Analyze the notes in the measure
    left_hand_notes = "G"
    right_hand_chord_tones = "G, B-natural, D"
    right_hand_non_chord_tone = "F#"
    
    print(f"Step 2: In measure {measure_number}, the left hand plays '{left_hand_notes}' octaves.")
    print(f"         The right hand plays chord tones '{right_hand_chord_tones}'.")
    print(f"         The '{right_hand_non_chord_tone}' is a brief, non-chord tone (an upper neighbor to G).")
    
    # Step 3: Determine the chord
    chord_notes = "G, B-natural, D"
    chord_name = "G major triad"
    print(f"Step 3: The essential notes forming the harmony are {chord_notes}, which create a {chord_name}.")
    
    # Step 4: Determine the Roman numeral
    scale_degree = 5
    function = "Dominant"
    roman_numeral = "V"
    
    print(f"Step 4: In the key of {key}, G is the {scale_degree}th scale degree, which is the {function}.")
    print(f"         A major triad built on the dominant is represented by the Roman numeral '{roman_numeral}'.")
    
    # Final Answer
    print("\nFinal Answer:")
    print(f"The correct Roman numeral for measure {measure_number} is {roman_numeral}.")

analyze_measure()