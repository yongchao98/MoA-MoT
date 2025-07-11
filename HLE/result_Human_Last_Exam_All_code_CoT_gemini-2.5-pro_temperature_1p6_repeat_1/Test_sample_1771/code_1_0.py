def trace_etymology():
    """
    Traces the etymology of a PIE verb through Germanic languages.
    """
    pie_form = "*seh₂gieti"
    meaning = "'(s)he’s giving a sign / seeking'"
    pie_root = "*sāg-" # The root, clarified from *seh₂g-

    # Step 1: Proto-Germanic
    # PIE *sāg- > PGmc *sōk- (Grimm's Law: g > k; Vowel shift: ā > ō)
    # PIE 3sg present *-eti > PGmc *-iþi / *-iþ
    # The verb is a class 1 weak verb (*sōkjaną). 3sg present is *sōkiþ.
    pgmc_form = "*sōkiþ"

    # Step 2: Old Norse
    # PGmc *sōkiþ > ON sækir
    # i-umlaut: *ō > œ/æ (caused by the *i of the stem)
    # Ending: *-iþ > -ir
    on_form = "sækir"

    # Step 3: Old English
    # PGmc *sōkiþ > OE sēceþ
    # i-umlaut: *ō > ē
    # Palatalization: *k > ċ (/tʃ/) before front vowel i
    # Ending: *-iþ > -eþ
    oe_form = "sēceþ"

    # Step 4: Old High German
    # PGmc *sōkiþ > OHG suohhit
    # Vowel shift: *ō > uo
    # High German Consonant Shift: *k > hh (/x/)
    # Ending: *-iþ > -it
    ohg_form = "suohhit"

    # Print the results in a clear, chained format
    print(f"The PIE form {pie_form}, with root {pie_root}, meaning {meaning}, develops as follows:")
    print("-" * 60)
    print(f"Proto-Indo-European: *seh₂gieti")
    print("    ↓")
    print(f"Proto-Germanic: {pgmc_form}")
    print("    ↓")
    print(f"    ├─ Old Norse: {on_form}")
    print(f"    ├─ Old English: {oe_form}")
    print(f"    └─ Old High German: {ohg_form}")

trace_etymology()