def analyze_measure_30():
    """
    Analyzes measure 30 of Beethoven's "Pathetique" Sonata, 1st movement,
    to determine the correct Roman numeral.
    """
    key = "C minor"
    measure = 30
    
    # Step 1: Identify the notes in the measure.
    left_hand_notes = ["G", "B-natural", "D"]
    right_hand_notes = ["G", "B-natural", "D"] # Harmonic tones
    chord_notes = sorted(list(set(left_hand_notes + right_hand_notes)))
    
    # Step 2: Identify the chord and its quality.
    chord_root = "G"
    chord_quality = "Major"
    chord_name = f"{chord_root} {chord_quality}"
    
    # Step 3: Relate the chord to the key of C minor.
    scale_degree = 5
    
    # Step 4: Determine the Roman numeral.
    # In music theory, a major chord built on the 5th degree (the dominant)
    # is represented by an uppercase Roman numeral 'V'.
    # This is standard practice in minor keys, using the raised 7th (leading tone)
    # from the harmonic minor scale.
    roman_numeral = "V"
    
    print(f"Analysis of Beethoven 'Pathetique' Sonata, 1st mvt., measure {measure}:")
    print("-" * 60)
    print(f"1. Key Signature: {key} (3 flats: B-flat, E-flat, A-flat).")
    print(f"2. Notes in measure {measure}: The primary notes are {', '.join(chord_notes)}.")
    print(f"3. Chord Identified: These notes form a {chord_name} triad.")
    print(f"4. Harmonic Function: In the key of {key}, the note {chord_root} is the {scale_degree}th scale degree (the Dominant).")
    print(f"5. Roman Numeral: A major triad built on the dominant is represented by the uppercase Roman numeral 'V'.")
    print("\nTherefore, the correct Roman numeral is:")
    print(roman_numeral)

analyze_measure_30()