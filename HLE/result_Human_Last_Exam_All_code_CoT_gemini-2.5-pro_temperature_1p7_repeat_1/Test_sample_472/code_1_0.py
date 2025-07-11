def scan_latin_line():
    """
    Scans a line of Latin, providing the metrical quantity (L/S) for each syllable.
    This line is from Plautus and follows the rules of dramatic meter, which can be irregular.
    The output shows the syllable-by-syllable analysis.
    """
    line = "et tibi bene esse soli quom sibi sit male"
    
    # The analysis of each word/elided pair based on classical rules:
    # 1. et: Long by position (followed by 't' in 'tibi').
    # 2. tibi: 'ti' is short, 'bi' is short.
    # 3. bene esse: elides to 'ben'esse'. Syllables are 'be' (short), 'nes' (long by position before 'ss'), 'se' (short).
    # 4. soli: Dative sōlī. 'sō' is long by nature, 'lī' is long by nature.
    # 5. quom: Long by position (followed by 's' in 'sibi').
    # 6. sibi: 'si' is short, 'bi' is short.
    # 7. sit: Long by position (followed by 'm' in 'male').
    # 8. male: 'ma' is short. The final syllable 'le' is anceps (treated as long at the end of a line).

    words = ["et", "tibi", "ben'esse", "soli", "quom", "sibi", "sit", "male"]
    scansion = ["L", "S S", "S L S", "L L", "L", "S S", "L", "S L"]

    print(f"Original line: \"{line}\"")
    print("\nScansion (syllable by syllable, grouped by word):\n")

    full_scansion_line = []
    full_word_line = []

    for i in range(len(words)):
        word = words[i]
        scan = scansion[i]
        print(f"{word.ljust(10)} -> {scan}")
        # Build a single line representation for the final output
        syllables = scan.split()
        full_scansion_line.extend(syllables)
    
    # To satisfy the format of having a final equation with numbers (syllables),
    # we will reconstruct the line with spaces between logical metrical units (feet).
    # Note: As this is not a regular hexameter, the foot division is an interpretation
    # of the complex comedic meter (Iambic Senarius with many irregularities).
    # Based on the analysis L S S S L S L L L S S L S L, a possible (but difficult) scansion is:
    # et ti-bi | ben' es- | se so- | li quom | si-bi | sit ma-le
    # L S S (Dactyl) | S L (Iamb) | S L (Iamb) | L L (Spondee) | S S (Pyrrhic) | L S L (Bacchius)
    # This shows the irregularity. We will insert spaces to reflect this interpretation.

    final_scansion_string = " ".join([
        "L S S", # et tibi
        "S L",   # ben' es-
        "S L",   # -se so-
        "L L",   # -li quom
        "S S",   # sibi
        "L S L"  # sit male
    ])

    print("\n---")
    print("Final scansion pattern with interpretive foot divisions:")
    # The final print has to output each part of the final "equation".
    print("et tibi bene esse soli quom sibi sit male")
    print(final_scansion_string)


scan_latin_line()