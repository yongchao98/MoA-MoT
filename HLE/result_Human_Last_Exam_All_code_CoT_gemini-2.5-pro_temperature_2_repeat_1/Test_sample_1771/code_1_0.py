import textwrap

def trace_etymology():
    """
    Traces a PIE word through its Germanic descendants based on linguistic sound laws.
    """
    # 1. Start with the given Proto-Indo-European (PIE) form.
    pie_form = "*seh₂gieti"
    pie_root = "*seh₂g-"
    pie_meaning = "'(s)he’s giving a sign'"
    print(f"Starting with PIE: {pie_form} {pie_meaning}")
    print("-" * 30)

    # 2. PIE to Proto-Germanic (*sōkijaną)
    pgmc_explanation = """
    From PIE to Proto-Germanic (PGmc), we see these key changes:
    1. Laryngeal Coloring: The laryngeal 'h₂' colors the preceding 'e' to 'a'. Root becomes '*sag-'.
    2. Vowel Shift: The resulting long vowel PIE *ā (from *eh₂) becomes PGmc *ō. Root is now '*sōk-'.
    3. Grimm's Law: The voiced stop *g becomes a voiceless stop *k.
    4. Suffix: The PIE suffix *-yeti becomes the PGmc Class 1 weak verb infinitive suffix *-janą.
    """
    pgmc_form = "*sōkijaną"
    pgmc_meaning = "'to seek'"
    print("Stage 1: Proto-Germanic")
    print(textwrap.dedent(pgmc_explanation).strip())
    print(f"Result in Proto-Germanic: *{pie_root} > {pgmc_form} ({pgmc_meaning})")
    print("-" * 30)

    # 3. Proto-Germanic to Old Norse (sœkja)
    on_explanation = """
    From PGmc to Old Norse (ON), the main change is i-umlaut:
    1. i-Umlaut: The 'j' in the suffix *-janą palatalizes (fronts) the preceding back vowel *ō to œ.
    2. Suffix Change: The PGmc infinitive ending evolves into the Old Norse '-a', and the 'j' is retained.
    """
    on_form = "sœkja"
    on_meaning = "'to seek, visit, prosecute'"
    print("Stage 2: Old Norse")
    print(textwrap.dedent(on_explanation).strip())
    print(f"Result in Old Norse: {pgmc_form} > {on_form} ({on_meaning})")
    print("-" * 30)

    # 4. Proto-Germanic to Old English (sēcan)
    oe_explanation = """
    From PGmc to Old English (OE), two crucial West Germanic changes occur:
    1. i-Umlaut: Like in Norse, the 'j' mutates the root vowel, but in OE, *ō becomes ē.
    2. Palatalization: The consonant *k is palatalized by the 'j' to ċ (pronounced /tʃ/, like 'ch' in 'church').
    3. The 'j' is then lost and the ending becomes '-an'.
    """
    oe_form = "sēcan"
    oe_meaning = "'to seek'"
    print("Stage 3: Old English")
    print(textwrap.dedent(oe_explanation).strip())
    print(f"Result in Old English: {pgmc_form} > {oe_form} ({oe_meaning})")
    print("-" * 30)

    # 5. Proto-Germanic to Old High German (suohhen)
    ohg_explanation = """
    From PGmc to Old High German (OHG), the changes are distinct:
    1. High German Consonant Shift: The voiceless stop *k shifts to a voiceless fricative, written 'hh' /x/.
    2. Diphthongization: The long vowel *ō is broken into the diphthong 'uo'.
    3. The suffix *-janą becomes '-en'. The 'j' causes an i-umlaut on 'uo' -> 'üe', but this was not always written.
    """
    ohg_form = "suohhen"
    ohg_meaning = "'to seek, to visit'"
    print("Stage 4: Old High German")
    print(textwrap.dedent(ohg_explanation).strip())
    print(f"Result in Old High German: {pgmc_form} > {ohg_form} ({ohg_meaning})")
    print("-" * 30)
    
    # Final Summary
    print("Summary of Etymological Development:")
    print(f"PIE *seh₂gieti gives us:")
    print(f"Proto-Germanic: {pgmc_form}")
    print(f"Old Norse: {on_form}")
    print(f"Old English: {oe_form}")
    print(f"Old High German: {ohg_form}")

if __name__ == "__main__":
    trace_etymology()