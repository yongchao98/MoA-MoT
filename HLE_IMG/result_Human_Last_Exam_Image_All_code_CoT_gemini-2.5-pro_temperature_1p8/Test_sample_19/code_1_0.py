def music_theory_analysis():
    """
    Provides a detailed music theory analysis for the provided excerpt of
    Beethoven's Moonlight Sonata.
    """

    # --- Part 1: Key of Modulation ---
    part1_title = "1. To what key the music modulates in measure 11 and the beginning of measure 12?"
    part1_answer = (
        "In measure 11, the chords are G dominant 7th (G-B-D-F) alternating with C major (C-E-G). "
        "This progression leads to and solidifies a C major chord on the first beat of measure 12. "
        "The G7 chord acts as the dominant (V7) resolving to its tonic (I). "
        "Therefore, the music modulates to the key of C major."
    )

    # --- Part 2: Justification for Modulation ---
    part2_title = "\n2. Given the 4 sharps key environment, what is the connection and/or justification for this particular modulation?"
    part2_answer = (
        "The home key of the piece is C# minor. A common modulation is to the relative major key, which is E major (the III chord), a key area established in measures 5-9.\n"
        "From this context of E major, the modulation to C major is a modulation to the flat-submediant (bVI). "
        "The diatonic submediant (vi) in E major is C# minor. Beethoven chromatically alters this by lowering the root C# to C natural, "
        "creating the C major chord. This movement to a chromatic mediant (bVI) is a powerful and colorful harmonic device characteristic of Beethoven and the Romantic era."
    )

    # --- Part 3: Roman Numeral Marking ---
    part3_title = "\n3. What would be the complete and accurate Roman numeral marking for the first beat of measure 11?"
    part3_answer = (
        "The chord on the first beat of measure 11 is a G dominant 7th chord (G7). Its function is to prepare and lead into the new key of C major. "
        "As established in Part 2, the C major chord functions as the flat-submediant (bVI) of the local key of E major. "
        "Therefore, the G7 chord is a secondary dominant. The complete and accurate Roman numeral marking is V7/bVI "
        "(read as 'five-seven of flat-six')."
    )
    
    final_numeral_equation = "V7 / bVI"
    
    print(part1_title)
    print(part1_answer)
    print(part2_title)
    print(part2_answer)
    print(part3_title)
    print(part3_answer)
    
    # Final combined answer in the requested format
    # The instructions asked to "output each number in the final equation"
    # I am interpreting this to mean spelling out the Roman numeral analysis
    full_answer_summary = f"""
1. The modulation is to C major.
2. The modulation is justified as a move to the flat-submediant (bVI) of the relative major key (E major).
3. The Roman numeral is V7/bVI. Broken down, this is V (5), 7, b (flat), VI (6)."""


if __name__ == '__main__':
    music_theory_analysis()
    final_answer = "1. C major.\n2. It is a modulation to the flat-submediant (bVI) of the relative major (E major).\n3. V7/bVI"
    print(f"\n<<<Combined Answer:\n{final_answer}>>>")
