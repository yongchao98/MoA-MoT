def analyze_beethoven_chord():
    """
    This function explains the harmonic analysis of measure 30
    of Beethoven's Pathetique Sonata, 1st movement, and prints the result.
    """
    key = "C minor"
    measure_number = 30
    
    # Step 1: Identify the notes in the measure
    notes_lh = ["G2", "G3"]
    notes_rh = ["G4", "B-natural4", "D5", "F-sharp5"]
    all_notes = "G, B-natural, D, F-sharp"
    
    # Step 2: Determine the root and function of the chord
    root_note = "G"
    scale_degree = 5
    function = "Dominant"
    numeral_base = "V"
    
    # Step 3: Determine the quality of the chord
    # The triad G-B(natural)-D is a major triad.
    # In a minor key, the dominant is typically major, so we use an uppercase V.
    quality = "Major"
    final_numeral = "V"

    # Step 4: Print the explanation and the final answer
    print(f"Analysis of Beethoven's Pathetique Sonata, 1st Mvt, Measure {measure_number}:")
    print(f"1. Key Signature: {key}")
    print(f"2. Notes in measure {measure_number}: {all_notes}")
    print(f"3. Bass Note/Root: {root_note}")
    print(f"4. The root '{root_note}' is the {scale_degree}th scale degree of {key}, making it the {function} chord.")
    print(f"5. The chord's triad (G-B-D) is {quality}, so the Roman numeral is uppercase.")
    print("6. Although the full sonority is a Major 7th chord (Vmaj7) due to the F#, its primary harmonic function is simply that of the Dominant.")
    print("\nTherefore, the correct Roman numeral is:")
    print(final_numeral)

analyze_beethoven_chord()