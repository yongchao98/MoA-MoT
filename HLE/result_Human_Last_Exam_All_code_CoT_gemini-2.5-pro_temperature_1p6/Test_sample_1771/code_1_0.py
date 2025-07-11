def trace_etymology():
    """
    This function traces the PIE root *seh₂g- through several Germanic languages
    by explaining the historical sound changes at each step.
    """
    # Define the forms for easy reference and printing.
    pie_verb = "*seh₂gieti"
    pie_root = "*sēh₂g-"
    pgmc_form = "*sōkijaną"
    on_form = "sœkja"
    oe_form = "sēcan"
    ohg_form = "suohhen"

    print(f"Tracing the etymology of PIE '{pie_verb}':\n")

    # 1. From Proto-Indo-European to Proto-Germanic
    print("1. From PIE to Proto-Germanic (PGmc):")
    print(f"  - The starting point is the PIE root '{pie_root}', inferred from cognates with a long vowel like Latin 'sāgīre'.")
    print("  - The PIE laryngeal '*h₂' colors the preceding vowel to '*a' and is then lost. PIE '*ē' thus becomes '*ā'.")
    print("  - PIE long '*ā' becomes PGmc long '*ō'.")
    print("  - By Grimm's Law, the PIE voiced stop '*g' becomes the PGmc voiceless stop '*k'.")
    print("  - The root evolves: *sēh₂g-  -> *sāg- -> *sōg- -> *sōk-.")
    print("  - The verb joins the first class of weak verbs, whose infinitive ends in *-janą.")
    print(f"  - Final Proto-Germanic form: '{pgmc_form}' (to seek).\n")

    # 2. From Proto-Germanic to Old Norse
    print(f"2. From PGmc '{pgmc_form}' to Old Norse (ON):")
    print("  - The 'j' in the suffix causes i-Umlaut, fronting the root vowel 'ō' to 'œ'.")
    print("  - The ending becomes '-a'.")
    print(f"  - Final Old Norse form: '{on_form}'.\n")

    # 3. From Proto-Germanic to Old English
    print(f"3. From PGmc '{pgmc_form}' to Old English (OE):")
    print("  - The 'j' in the suffix causes i-Umlaut, fronting the root vowel 'ō' to 'ē'.")
    print("  - The 'j' is then lost after the long vowel, and the ending becomes '-an'.")
    print(f"  - Final Old English form: '{oe_form}'.\n")
    
    # 4. From Proto-Germanic to Old High German
    print(f"4. From PGmc '{pgmc_form}' to Old High German (OHG):")
    print("  - The High German Consonant Shift changes the voiceless stop '*k' to the long fricative 'hh' [xx].")
    print("  - PGmc 'ō' becomes the diphthong 'uo' in Old High German.")
    print("  - The standard written form often doesn't show i-Umlaut ('uo' -> 'üe'). The ending becomes '-en'.")
    print(f"  - Final Old High German form: '{ohg_form}'.\n")

    # Final summary with all derived forms
    print("--- Summary ---")
    print(f"PIE '{pie_verb}' gives us the following infinitive forms:")
    print(f"1. Proto-Germanic: {pgmc_form}")
    print(f"2. Old Norse: {on_form}")
    print(f"3. Old English: {oe_form}")
    print(f"4. Old High German: {ohg_form}")


if __name__ == '__main__':
    trace_etymology()