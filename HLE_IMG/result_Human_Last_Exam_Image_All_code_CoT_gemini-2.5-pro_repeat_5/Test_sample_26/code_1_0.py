def analyze_beethoven_measure():
    """
    Determines the Roman numeral for the chord in measure 30 of
    Beethoven's "Pathetique" Sonata, first movement.
    """
    # Step 1: Identify the key signature.
    # The key signature has 3 flats (Bb, Eb, Ab), and the piece's tonic is C.
    key = "C minor"
    measure = 30
    print(f"Step 1: The musical key for measure {measure} is {key}.")

    # Step 2: Identify the notes in the measure.
    # The right hand plays an arpeggio of G, B-natural, D, F-sharp.
    # The left hand plays a G octave.
    notes = ["G", "B-natural", "D", "F-sharp"]
    print(f"Step 2: The notes found in measure {measure} are: {', '.join(notes)}.")

    # Step 3: Determine the root and quality of the chord.
    # The notes G, B, D, and F# form a G major seventh chord.
    root = "G"
    quality = "Major Seventh"
    chord_name = "Gmaj7"
    print(f"Step 3: The chord is a {root} {quality} ({chord_name}).")

    # Step 4: Determine the function of the chord within the key.
    # In C minor, G is the 5th scale degree, which is the dominant.
    scale_degree = 5
    function = "Dominant"
    print(f"Step 4: The chord's root, G, is the {function} (scale degree {scale_degree}) of {key}.")

    # Step 5: Construct the Roman numeral.
    # A major chord on the dominant (5th degree) is written with an uppercase 'V'.
    # A major seventh is added to the numeral as 'maj7'.
    roman_numeral_base = "V"
    roman_numeral_quality = "maj7"
    final_numeral = roman_numeral_base + roman_numeral_quality
    print(f"Step 5: The Roman numeral for a {quality} chord on the degree {scale_degree} is '{final_numeral}'.")
    print("\n--- FINAL ANSWER ---")
    print(f"The correct Roman numeral is: {final_numeral}")

analyze_beethoven_measure()