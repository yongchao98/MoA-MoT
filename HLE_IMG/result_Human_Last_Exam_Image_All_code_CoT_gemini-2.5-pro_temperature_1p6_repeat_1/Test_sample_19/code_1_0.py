def analyze_moonlight_sonata_modulation():
    """
    Provides a step-by-step music theory analysis for the specified
    passage in Beethoven's Moonlight Sonata, 1st Movement.
    """
    
    print("--- Moonlight Sonata Modulation Analysis (Measures 11-12) ---")

    # Part 1: The New Key
    part1_title = "\n1. To what key does the music modulate in measure 11 and the beginning of measure 12?"
    part1_answer = (
        "The music modulates to the key of B Major.\n"
        "This is established by the clear B Major chord (notes B-D#-F#) on the first beat of measure 12. "
        "This new tonic is prepared by its dominant chord, F# Major (notes F#-A#-C#), throughout the second half of measure 11."
    )
    print(part1_title)
    print(part1_answer)

    # Part 2: The Justification
    part2_title = "\n2. Given the 4-sharp key environment, what is the justification for this modulation?"
    part2_answer = (
        "The original key of the piece is C# minor. Its relative major key is E Major.\n"
        "The new key, B Major, is the Dominant (the V chord) of E Major.\n"
        "Therefore, the modulation is from the tonic minor (i, C# minor) to the 'Dominant of the Relative Major' (V of III). "
        "This is a common and sophisticated harmonic device used to create contrast and build tension while remaining logically connected to the home key."
    )
    print(part2_title)
    print(part2_answer)
    
    # Part 3: The Roman Numeral
    part3_title = "\n3. What is the complete and accurate Roman numeral marking for the first beat of measure 11?"
    part3_answer = (
        "The chord on the first beat of measure 11 contains the notes G#, B, and D-natural over a B in the bass.\n"
        "This forms a G# diminished triad (G#-B-D) in first inversion.\n\n"
        "Analyzed in the home key of C# minor:\n"
        "- A G# diminished triad is the leading-tone chord (vii°) to the key of A minor/major.\n"
        "- The key of 'A' is the submediant (the VI chord) of C# minor.\n"
        "- Thus, this chord functions as a secondary leading-tone chord tonicizing the submediant area.\n"
        "- Since the chord is in first inversion (with the '6' indicating the interval of a sixth above the bass), the final and complete Roman numeral is vii°⁶/VI."
    )
    print(part3_title)
    print(part3_answer)

    # Final summary of the Roman Numeral parts
    final_numeral = "vii°⁶/VI"
    print(f"\nFinal Roman Numeral Marking: {final_numeral}")
    print("Broken down: vii° (leading-tone diminished) ⁶ (in first inversion) /VI (of the submediant)")


# Execute the analysis function
analyze_moonlight_sonata_modulation()