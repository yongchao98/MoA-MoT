def trace_etymology():
    """
    This script traces a PIE root through its development into various Germanic languages.
    """

    # 1. The Proto-Indo-European (PIE) starting point.
    pie_form = "*seh₂gieti"
    pie_meaning = "(s)he’s giving a sign"

    # 2. Proto-Germanic (PGmc)
    # The attested form *sōkijaną 'to seek' develops from the o-grade of the PIE root (*soh₂g-).
    # Grimm's Law: PIE *g -> PGmc *k.
    # Vowel Change: PIE *o + *h₂ -> PGmc *ō.
    pgmc_form = "*sōkijaną"

    # 3. Old Norse (ON)
    # i-umlaut: The 'j' in *sōkijaną causes the root vowel *ō to become œ.
    # The infinitive ending *-aną becomes -a.
    on_form = "sœkja"

    # 4. Old English (OE)
    # i-umlaut: The 'j' causes *ō > œ̄ > ē. The 'j' is then lost.
    # The infinitive ending *-aną becomes -an.
    oe_form = "sēcan"

    # 5. Old High German (OHG)
    # High German Consonant Shift: PGmc *k -> OHG hh (or ch).
    # Vowel Breaking: PGmc *ō -> OHG uo.
    # The infinitive ending *-aną becomes -an.
    ohg_form = "suohhan"

    # 6. Print the full etymological chain.
    print(f"PIE {pie_form} > Proto-Germanic {pgmc_form} > Old Norse {on_form} > Old English {oe_form} > Old High German {ohg_form}")

trace_etymology()