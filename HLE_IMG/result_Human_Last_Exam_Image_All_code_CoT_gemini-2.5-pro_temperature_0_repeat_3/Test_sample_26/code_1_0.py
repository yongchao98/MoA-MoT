def analyze_measure_30():
    """
    Analyzes measure 30 of Beethoven's "Pathetique" Sonata, 1st mvt,
    to determine the correct Roman numeral.
    """
    key = "C minor"
    
    # Step 1: Identify the notes in measure 30
    # Treble clef (right hand): G, B-natural, D, F-sharp
    # Bass clef (left hand): G, D
    notes_in_measure = ["G", "B-natural", "D", "F-sharp"]
    
    # Step 2: Determine the chord
    # The notes G, B-natural, and D form a G major triad.
    # The F-sharp is the major seventh above the root G.
    # The full chord is a G major seventh (Gmaj7).
    chord_root = "G"
    chord_quality = "Major"
    
    # Step 3: Analyze the chord's function in the key of C minor
    # The scale degrees of C minor are:
    # 1: C (tonic)
    # 2: D
    # 3: Eb
    # 4: F
    # 5: G (dominant)
    # 6: Ab
    # 7: Bb (or B-natural for harmonic/melodic minor)
    
    # The root of the chord is G, which is the 5th degree of C.
    # Therefore, the chord has a dominant function.
    function = "Dominant"
    scale_degree = 5
    
    # Step 4: Assign the Roman Numeral
    # A dominant chord on the 5th scale degree is represented by the Roman numeral V.
    # In a minor key, the dominant chord is typically major (using the raised 7th, or leading tone),
    # which is represented by an uppercase V.
    # The F-sharp is a chromatic alteration, but the fundamental harmony is the dominant (V).
    roman_numeral = "V"
    
    print("Analysis of Beethoven's 'Pathetique' Sonata, Op. 13, 1st Mvt, Measure 30:")
    print("-" * 70)
    print(f"Key Signature: {key}")
    print(f"Notes in Measure: {', '.join(notes_in_measure)}")
    print(f"Chord: {chord_root} {chord_quality} (with a chromatic F#)")
    print(f"Function: The chord is built on the {scale_degree}th scale degree, giving it a {function} function.")
    print("\nThe chord G-B(natural)-D is the standard dominant triad in C minor.")
    print("It strongly leads back to the tonic chord (C minor) in the following measure.")
    print("The F-sharp is a chromatic coloration but does not change the chord's fundamental dominant function.")
    print("\nTherefore, the correct Roman numeral is:")
    print(f"The final Roman numeral is: {roman_numeral}")

analyze_measure_30()