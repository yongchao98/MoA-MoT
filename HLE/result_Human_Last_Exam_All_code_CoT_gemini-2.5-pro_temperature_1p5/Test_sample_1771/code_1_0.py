def solve_linguistic_derivation():
    """
    This function traces the derivation of a Proto-Indo-European root
    through Proto-Germanic into its daughter languages and prints the result.
    """

    pie_form = "*seh₂gieti"
    pie_meaning = "'(s)he gives a sign, tracks'"
    
    pgmc_form = "*sōkijaną"
    on_form = "sœkja"
    oe_form = "sēcan"
    ohg_form = "suohhen"
    
    print("This shows the derivation of a PIE root into the early Germanic languages.")
    print("-" * 70)
    print(f"1. Proto-Indo-European (PIE): {pie_form}")
    print(f"   Meaning: {pie_meaning}")
    print(f"   Root: *seh₂g-\n")

    print("2. Sound Changes to Proto-Germanic (PGmc):")
    print("   - Grimm's Law: Voiced stop *g > voiceless stop *k.")
    print("   - Laryngeal Effect: *eh₂ > *ā.")
    print("   - PGmc Vowel Shift: *ā > *ō.")
    print("   - Suffix Change: PIE *-ie- > PGmc *-ja- (Class 1 weak verb).")
    print(f"   - Resulting PGmc Infinitive: {pgmc_form} ('to seek')\n")
    
    print("3. Development into Daughter Languages:\n")

    # Old Norse
    print(f"   - Proto-Germanic: {pgmc_form}")
    print(f"     -> Old Norse: {on_form}")
    print(f"        (Sound changes: i-umlaut of *ō > œ)\n")

    # Old English
    print(f"   - Proto-Germanic: {pgmc_form}")
    print(f"     -> Old English: {oe_form}")
    print(f"        (Sound changes: i-umlaut of *ō > ē; palatalization of *k > ċ /tʃ/)\n")

    # Old High German
    print(f"   - Proto-Germanic: {pgmc_form}")
    print(f"     -> Old High German: {ohg_form}")
    print(f"        (Sound changes: PGmc *ō > uo; High German Shift of *k > hh /x/)\n")

    print("-" * 70)
    print("Final Summary:")
    print(f"Proto-Germanic: {pgmc_form}")
    print(f"Old Norse: {on_form}")
    print(f"Old English: {oe_form}")
    print(f"Old High German: {ohg_form}")

solve_linguistic_derivation()