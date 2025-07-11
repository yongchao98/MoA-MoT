def analyze_measure_30():
    """
    Analyzes the harmony of measure 30 from Beethoven's Pathetique Sonata, Mvt. 1,
    and provides the correct Roman numeral.
    """
    piece = "Beethoven's 'Pathetique' Sonata, Op. 13, 1st Movement"
    measure = 30
    home_key = "C minor"

    print(f"Analyzing Measure {measure} of {piece}")
    print("-" * 50)
    print(f"1. Home Key Signature: {home_key} (B♭, E♭, A♭)")

    # Notes in measure 30
    rh_notes = ["F♯", "A", "C", "E♭"]
    lh_notes = ["D"]
    all_notes = sorted(list(set(lh_notes + rh_notes)), key="CDEFGAB".index)

    print("\n2. Notes in Measure 30:")
    print(f"   - Right Hand (Treble Clef): {', '.join(rh_notes)}")
    print(f"   - Left Hand (Bass Clef): {', '.join(lh_notes)}")
    print(f"   - Combined Notes: {', '.join(all_notes)}")

    # Chord Identification
    root = "D"
    quality = "Dominant 7th with a minor 9th (D7♭9)"
    chord_notes_explained = "D (root), F♯ (major third), A (perfect fifth), C (minor seventh), E♭ (minor ninth)"

    print("\n3. Chord Identification:")
    print(f"   - Root: {root}")
    print(f"   - Chord Quality: {quality}")
    print(f"   - Spelled out: {chord_notes_explained}")

    # Harmonic Function
    resolves_to = "G minor (in measure 31)"
    function_in_c_minor = "The dominant chord of G minor."
    g_minor_in_c_minor = "v (the minor dominant)"
    secondary_dominant_name = "V7/v ('five-seven of five')"

    print("\n4. Harmonic Function:")
    print(f"   - The chord in measure 30 resolves to a {resolves_to}.")
    print(f"   - A D7 chord is the dominant of G. In our home key of {home_key}, G minor is the 'v' chord.")
    print("   - Therefore, the D7(♭9) chord is a secondary dominant.")

    # Final Roman Numeral
    final_roman_numeral = "V7/v"

    print("\n5. Conclusion:")
    print("   The Roman numeral describes both the chord's quality (Dominant 7th) and its function (secondary dominant).")
    print(f"\nThe correct Roman numeral for measure {measure} is: {final_roman_numeral}")

analyze_measure_30()