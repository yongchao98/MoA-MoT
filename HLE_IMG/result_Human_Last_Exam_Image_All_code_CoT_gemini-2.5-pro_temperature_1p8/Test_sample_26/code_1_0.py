def analyze_pathetique_measure():
    """
    Analyzes measure 30 of Beethoven's Pathetique Sonata (Op. 13, Mvt. 1)
    to determine its Roman numeral.
    """
    key = "C minor"
    measure_number = 30
    
    print("Step 1: Identify the Key Signature")
    print(f"The key signature is three flats (Bb, Eb, Ab). The piece is in {key}.")
    print("-" * 30)

    print(f"Step 2: Identify the Notes in Measure {measure_number}")
    notes_in_measure = ["G#", "B", "D", "F"]
    print("The notes played throughout measure 30 are G#, B, D, and F.")
    print("-" * 30)

    print("Step 3: Identify the Chord")
    chord_name = "G-sharp fully diminished seventh (G#°7)"
    print(f"These notes form a {chord_name}.")
    print("A diminished seventh chord is built with a root, a minor third, a diminished fifth, and a diminished seventh.")
    print("-" * 30)

    print("Step 4: Analyze the Harmonic Function in C minor")
    leading_tone = "B"
    print(f"The key is C minor. The leading tone (the 7th scale degree) of C minor is {leading_tone}.")
    print("A leading-tone chord is built on this 7th degree and resolves to the tonic (C minor).")
    
    leading_tone_chord_notes = "B - D - F - Ab"
    print(f"The leading-tone diminished seventh chord in C minor (vii°7) is built from the notes: {leading_tone_chord_notes}.")
    
    print("\nThe notes in measure 30 are G# - B - D - F.")
    print("Notice that G# is the same pitch as Ab (they are enharmonic equivalents).")
    print("Therefore, the G#°7 chord in measure 30 is functioning as the leading-tone seventh chord of C minor.")
    print("Beethoven spelled the chord with a G# for smoother voice leading into the next measure.")
    print("-" * 30)

    print("Step 5: State the Final Roman Numeral")
    final_numeral = "vii°7"
    print("The Roman numeral for the leading-tone fully diminished seventh chord in a minor key is vii°7.")
    print(f"\nThe correct Roman numeral for measure {measure_number} is: {final_numeral}")

analyze_pathetique_measure()