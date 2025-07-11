def scan_hexameter_line():
    """
    Scans the Latin hexameter line: "verum ubi equi atque hominis casu convenit imago."
    This function breaks down the scansion foot by foot based on the standard reading.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    scansion_pattern = "D D S D D S"
    
    print(f"Latin Hexameter Line: {line}")
    print(f"Standard Scansion: {scansion_pattern}\n")
    print("--- Foot-by-Foot Analysis ---")

    analysis = [
        "Foot 1 (Dactyl): 'verum ubi'",
        "   Syllables: vē - rum ubi -> vē - r'u - bi",
        "   Analysis: LONG - SHORT - SHORT.",
        "   - 'vē-' is long by nature (from vērus).",
        "   - '-rum' elides with the 'u' of 'ubi'.",
        "   - 'u-' and 'bi-' are both short.",

        "Foot 2 (Dactyl): 'equi at-'",
        "   Syllables: ē - qui - at",
        "   Analysis: LONG - SHORT - SHORT.",
        "   - 'ē-' is read as long here by poetic license, a common but difficult feature.",
        "   - 'qui' is a short syllable.",
        "   - 'at-' from 'atque' serves as the second short syllable before the elision with 'hominis'.",

        "Foot 3 (Spondee): '-que hominis'",
        "   Syllables: que ho - mi- -> qu'hō - nis",
        "   Analysis: LONG - LONG.",
        "   - 'que' elides with 'hominis', leaving 'qu'hō-', which is long by position (o + m + n).",
        "   - 'nis' is long by position (i + n + s of 'casu').",

        "Foot 4 (Dactyl): 'casu con-'",
        "   Syllables: cā - su - con",
        "   Analysis: LONG - SHORT - SHORT.",
        "   - 'cā-' is long by nature (ablative case).",
        "   - 'su-' and 'con-' (from convenit) are treated as the two short syllables.",
        
        "Foot 5 (Dactyl): '-venit i-'",
        "   Syllables: vē - nit - i",
        "   Analysis: LONG - SHORT - SHORT.",
        "   - 'vē-' is long by nature (from venio, perfect tense).",
        "   - 'nit' and the first syllable 'i-' of 'imago' are short.",

        "Foot 6 (Spondee): 'mago.'",
        "   Syllables: mā - gō",
        "   Analysis: LONG - LONG.",
        "   - 'mā-' and '-gō' are both long by nature."
    ]

    for line in analysis:
        print(line)

    final_result = "verum u | bi equi | atque ho | minis ca | su conve | nit imago"
    print("\n--- Final Result ---")
    print("Final Pattern: D D S D D S")
    

scan_hexameter_line()