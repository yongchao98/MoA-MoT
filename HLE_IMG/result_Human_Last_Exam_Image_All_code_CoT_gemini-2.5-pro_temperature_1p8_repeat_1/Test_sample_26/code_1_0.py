def analyze_beethoven_pathetique_m30():
    """
    Provides a step-by-step harmonic analysis of measure 30
    from Beethoven's "Pathetique" Sonata, 1st movement, and prints
    the resulting Roman numeral.
    """
    
    # Step 1: Define the musical context
    piece = "Beethoven, Piano Sonata No. 8 'Pathetique', Op. 13, Mvt. 1"
    measure = 30
    
    # Step 2: Identify the notes and the local key
    # The notes in measure 30 are G#, B, C#, E.
    # The chord resolves to A minor in measure 31, establishing a local key.
    local_key = "A minor"
    tonic = "A"
    
    # Step 3: Analyze the chord's root and quality
    # The bass note is G#, which is the leading tone (7th degree) of A minor.
    # The structural harmony is a G-sharp diminished triad (G#-B-D).
    chord_root_note = "G#"
    scale_degree_number = 7
    scale_degree_name = "Leading Tone"
    chord_quality = "diminished"
    
    # Step 4: Construct the Roman numeral
    # A diminished triad on the leading tone is denoted by 'vii°'.
    roman_numeral_symbol = "vii°"

    print(f"Harmonic Analysis of {piece}, Measure {measure}:")
    print("-" * 50)
    print(f"1. Local Key Context: The harmony points to '{local_key}'.")
    print(f"2. Chord Root: The bass note is {chord_root_note}.")
    print(f"3. Functional Analysis: {chord_root_note} is scale degree {scale_degree_number} (the {scale_degree_name}) in the key of {local_key}.")
    print(f"4. Chord Quality: The chord is a {chord_quality} triad.")
    print("-" * 50)
    print(f"The resulting Roman numeral for the chord in measure {measure} is: {roman_numeral_symbol}")

analyze_beethoven_pathetique_m30()