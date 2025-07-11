def analyze_moonlight_sonata_excerpt():
    """
    Provides a detailed music theory analysis of the specified excerpt from
    Beethoven's Moonlight Sonata and prints the findings.
    """
    
    print("--- Analysis of Beethoven's Moonlight Sonata, 1st Mvt., Measures 11-12 ---\n")

    # --- Part 1: Key of Modulation ---
    part1_title = "Part 1: The Key of Modulation"
    part1_analysis = [
        "1. The home key of the piece is C# minor. The 4-sharp key signature (F#, C#, G#, D#) can indicate E Major or C# minor, and the piece's harmony is clearly centered on C# minor.",
        "2. In measure 11, the harmony shifts. Chords like the A# diminished on beat 3 act as a leading-tone chord pointing towards B minor.",
        "3. However, this B minor functions as a dominant (V) to the chord that arrives on the downbeat of measure 12.",
        "4. Measure 12 begins with a clear, root-position E minor chord (Notes: E, G-natural, B). This chord is the new point of harmonic stability.",
        "Conclusion: The music modulates to the key of E minor."
    ]
    
    print(part1_title)
    print("-" * len(part1_title))
    for line in part1_analysis:
        print(line)
    
    print("\n" + "="*70 + "\n")
    
    # --- Part 2: Justification for Modulation ---
    part2_title = "Part 2: Connection and Justification"
    part2_analysis = [
        "1. The home key is C# minor. A standard modulation from a minor key is to its relative major, which shares the same key signature. The relative major of C# minor is E Major (the III chord).",
        "2. Beethoven modulates to E minor (the iii chord), not E Major. The destination key shares the same root (E) as the relative major, but it is in the minor mode.",
        "3. This technique is a form of 'mode mixture'. It allows for a modulation to a closely related key center (the mediant, or third scale degree) while preserving the dark, somber character of the piece.",
        "Conclusion: The modulation to the minor mediant (E minor) is justified as a way to create harmonic interest while maintaining the established mood."
    ]

    print(part2_title)
    print("-" * len(part2_title))
    for line in part2_analysis:
        print(line)
        
    print("\n" + "="*70 + "\n")

    # --- Part 3: Roman Numeral Analysis ---
    part3_title = "Part 3: Roman Numeral for Measure 11, Beat 1"
    part3_analysis = [
        "1. The notes on the first beat of measure 11 are: B in the bass, with G-natural and E-natural above it.",
        "2. These notes (E, G, B) form an E minor chord.",
        "3. We analyze this chord in the home key of C# minor:",
        "   - The root, E, is the third scale degree (mediant). Since the chord is minor, we use a lowercase Roman numeral: iii.",
        "   - The bass note, B, is the 5th of the E minor chord. A chord with its 5th in the bass is in second inversion.",
        "   - Second inversion is indicated by the figured bass symbols 6/4.",
        "Conclusion: The complete Roman numeral is iii⁶₄ (minor mediant in second inversion)."
    ]
    
    print(part3_title)
    print("-" * len(part3_title))
    for line in part3_analysis:
        print(line)
    
    print("\n--- Final Roman Numeral Components ---")
    numeral = "iii"
    upper_figure = 6
    lower_figure = 4
    print(f"Numeral: {numeral}")
    print(f"Upper Figure: {upper_figure}")
    print(f"Lower Figure: {lower_figure}")

# Execute the analysis function
analyze_moonlight_sonata_excerpt()