def solve_music_theory():
    """
    Analyzes the provided musical excerpt from Beethoven's Moonlight Sonata
    and answers the three-part question.
    """

    # Part 1: Key of Modulation
    # The key signature is 4 sharps (C# minor).
    # Measure 11 contains a B minor chord (B-D natural-F#).
    # Measure 12 begins with an E major chord (E-G#-B).
    # The progression from B minor to E major is a v-I cadence in the key of E major.
    key_modulation = "E major"

    # Part 2: Justification for the Modulation
    # The home key is C# minor. The new key is E major.
    # E major is the relative major of C# minor. They share the same key signature.
    justification = "The home key of the piece is C# minor. The new key, E major, is the relative major of C# minor. Modulating to the relative major is a very common and natural-sounding progression in tonal music because the two keys share the same notes and key signature (4 sharps)."

    # Part 3: Roman Numeral Analysis
    # The chord on the first beat of measure 11 is B minor.
    # This chord functions as the minor dominant (v) of the new key, E major.
    # The new key, E major, is the mediant (III) chord in the home key of C# minor.
    # Therefore, the chord is the minor dominant of the mediant.
    roman_numeral = "v/III"
    explanation = "The chord is a B minor chord, which is the minor dominant (v) of E major. Since E major is the mediant (III) of the home key (C# minor), the complete analysis is v/III."

    print("--- Music Theory Analysis ---")
    print(f"1. In measures 11-12, the music modulates to the key of: {key_modulation}")
    print("\n2. Connection and Justification:")
    print(f"   {justification}")
    print(f"\n3. Roman Numeral for the first beat of measure 11: {roman_numeral}")
    print(f"   Explanation: {explanation}")

solve_music_theory()