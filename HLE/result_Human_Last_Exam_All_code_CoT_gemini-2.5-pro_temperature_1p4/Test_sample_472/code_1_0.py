def scan_latin_line():
    """
    Scans the provided Latin line and prints the metrical analysis.
    L = Long syllable
    S = Short syllable
    """

    # The line is an iambic senarius with several substitutions.
    # 1. ben(e) esse -> Elision removes one syllable.
    # 2. Final syllable of 'male' is treated as long (syllaba anceps).
    # 3. Substitutions: Dactyl (LSS), Spondee (LL), Tribrach (SSS).

    line = [
        ("et", "L"), ("ti", "S"), ("bi", "S"),  # Foot 1: Dactyl (L S S)
        ("ben'", "S"), ("es", "L"),              # Foot 2: Iamb (S L)
        ("se", "S"), ("so", "L"),                # Foot 3: Iamb (S L)
        ("li", "L"), ("quom", "L"),              # Foot 4: Spondee (L L)
        ("si", "S"), ("bi", "S"), ("sit", "S"),  # Foot 5: Tribrach (S S S)
        ("ma", "S"), ("le", "L")                 # Foot 6: Iamb (S L)
    ]

    feet_boundaries = [3, 5, 7, 9, 12, 14]
    
    syllable_line = []
    scansion_line = []
    
    for i, (syllable, mark) in enumerate(line):
        syllable_line.append(syllable)
        scansion_line.append(mark)
        if (i + 1) in feet_boundaries and (i + 1) != len(line):
            syllable_line.append("|")
            scansion_line.append("|")

    # To ensure proper alignment, we format the output
    # Find the max width for each column
    col_widths = [max(len(s), len(m)) for s, m in zip(syllable_line, scansion_line)]

    # Build the formatted strings
    formatted_syllables = " ".join([syllable.center(col_widths[i]) for i, syllable in enumerate(syllable_line)])
    formatted_scansion = " ".join([mark.center(col_widths[i]) for i, mark in enumerate(scansion_line)])

    print("Scansion of the line:")
    print("\"et tibi bene esse soli quom sibi sit male\"")
    print("-" * 50)
    print(formatted_syllables)
    print(formatted_scansion)
    print("-" * 50)

scan_latin_line()