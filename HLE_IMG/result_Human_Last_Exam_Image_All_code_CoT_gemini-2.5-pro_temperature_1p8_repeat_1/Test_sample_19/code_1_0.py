def analyze_moonlight_sonata_passage():
    """
    Provides a detailed music theory analysis of the specified passage
    from Beethoven's Moonlight Sonata.
    """

    # Part 1: Key Modulation
    answer1_header = "1. To what key the music modulates in measure 11 and the beginning of measure 12?"
    answer1_text = """   The original key signature (4 sharps) and opening chords establish the key as C# minor.
   The chord at the beginning of measure 12 consists of the notes E, G#, and B.
   This is a clear E Major triad, acting as a point of arrival.
   Therefore, the music modulates to E Major."""

    # Part 2: Justification
    answer2_header = "\n2. Given the 4 sharps key environment, what is the connection and / or justification for this particular modulation?"
    answer2_text = """   E Major is the relative major of C# minor.
   These two keys share the exact same key signature (F#, C#, G#, D#), making modulation between them sound smooth and natural.
   Modulating from a minor key to its relative major is a fundamental technique in the common practice period, often used to introduce a contrasting, brighter-sounding theme or section."""

    # Part 3: Roman Numeral
    answer3_header = "\n3. Building on 1 and 2 answers - what would be the complete and accurate roman numeral marking for the first beat of measure 11?"
    answer3_text = """   The analysis must be done in the context of the new key, E Major.
   The notes on the first beat of measure 11 are:
   - Bass: B
   - Harmony/Melody: B, D#, F-double-sharp (Fx)

   These notes form a B augmented triad (B-D#-Fx). In the key of E Major, the V (dominant) chord is B Major (B-D#-F#).
   The chord in question is an altered dominant where the fifth is raised, creating immense tension that resolves to the new tonic (E) in the next measure.
   An augmented triad built on the dominant is represented by a 'V' with a plus sign.
   The complete roman numeral marking is: V+"""

    print(answer1_header)
    print(answer1_text)
    print(answer2_header)
    print(answer2_text)
    print(answer3_header)
    print(answer3_text)

# Execute the analysis function
analyze_moonlight_sonata_passage()