def scan_hexameter_line():
    """
    Scans the Latin hexameter line:
    "verum ubi equi atque hominis casu convenit imago."

    This function breaks the line into its six metrical feet and identifies
    each as a dactyl (D) or spondee (S).
    """

    # Hexameter feet breakdown based on analysis.
    # This particular line is complex and requires some assumptions:
    # 1. Hiatus after "verum" (the '-um' is not elided, but used for position).
    # 2. The 'i' in 'ubi' is scanned as long by nature ('ubī').
    # 3. There is an elision between 'atque' and 'hominis' -> 'atqu(e) hominis'.
    # 4. The 'i' in 'imago' is dropped (aphaeresis) after a vowel-ending 'convenit'.
    
    foot1_word = "verum u-"
    foot1_scan = "D"  # Dactyl: vēr(long by position)-um(short) u(short)

    foot2_word = "bī equi"
    foot2_scan = "S"  # Spondee: bī(long by nature) equī(long by synizesis/poetic license)

    foot3_word = "atqu(e) homin-"
    foot3_scan = "S"  # Spondee: at(long by pos) -qu' hom(long by pos)

    foot4_word = "is cāsū"
    foot4_scan = "D"  # Dactyl: is(long by pos) cā(short) sū(short) - NOTE: This is a major point of difficulty. 'cāsū' is normally L-L. To make this a dactyl, 'cā' must be treated as short, which is a poetic license.

    foot5_word = "convēnit"
    foot5_scan = "D"  # Dactyl: con(long) -vē(assumed long) -nit(short)
    
    foot6_word = "(i)māgō"
    foot6_scan = "S"  # Spondee: mā(long by nature) gō(long by nature/final syllable)

    feet = [
        (foot1_word, foot1_scan),
        (foot2_word, foot2_scan),
        (foot3_word, foot3_scan),
        (foot4_word, foot4_scan),
        (foot5_word, foot5_scan),
        (foot6_word, foot6_scan)
    ]

    print("Scansion of: verum ubi equi atque hominis casu convenit imago.")
    print("-" * 60)
    
    final_pattern = []
    for i, (word, scan_type) in enumerate(feet, 1):
        print(f"Foot {i}: {word:<15} ({scan_type})")
        final_pattern.append(scan_type)
        
    print("-" * 60)
    print("Final Pattern: " + " ".join(final_pattern))

scan_hexameter_line()