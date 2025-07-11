def trace_etymology():
    """
    Traces the PIE word *seh₂gieti to Proto-Germanic, Old Norse, Old English, and Old High German.
    """
    pie_form = "*seh₂gieti"
    print(f"Starting with the Proto-Indo-European (PIE) form: {pie_form}\n")

    # --- Step 1: Proto-Germanic ---
    print("--- 1. Derivation to Proto-Germanic (PGmc) ---")
    print("The PIE form is a 3rd person singular yod-present verb meaning 's/he seeks/perceives'.")
    print(f"PIE: {pie_form}")
    print("  a) The laryngeal *h₂ colors the preceding *e to *a, and then PIE *ā becomes PGmc *ō.")
    print("     Root *seh₂g- -> *sāg- -> PGmc *sōk-")
    print("  b) Grimm's Law affects the consonants: *g > *k and the verbal ending *-ti > *-þi.")
    print("  c) The PIE yod-present (*-ie-) becomes the Germanic Class 1 weak verb stem (infinitive *-janą, 3sg present stem *-i-).")
    print("  d) Assembling the 3rd person singular present tense form:")
    pgmc_form = "*sōkiþ"
    print(f"Resulting Proto-Germanic form: {pgmc_form}\n")

    # --- Step 2: Old Norse ---
    print("--- 2. Derivation to Old Norse (ON) ---")
    print(f"Starting from PGmc: {pgmc_form}")
    print("  a) I-Umlaut: The *i in the stem causes the root vowel *ō to be fronted to 'œ'.")
    print("     *sōk- > sœk-")
    print("  b) The 3rd person singular ending *-iþ evolves into '-ir' in Old Norse.")
    print("     -iþ > -ir")
    on_form = "sœkir"
    print(f"Resulting Old Norse form: {on_form}\n")

    # --- Step 3: Old English ---
    print("--- 3. Derivation to Old English (OE) ---")
    print(f"Starting from PGmc: {pgmc_form}")
    print("  a) I-Umlaut: The *i in the stem causes the root vowel *ō to be fronted to 'ē'.")
    print("     *sōk- > sēk-")
    print("  b) Palatalization: The consonant *k is palatalized to 'ċ' (pronounced /tʃ/) before the front vowel 'e'.")
    print("     sēk- > sēċ-")
    print("  c) The ending *-iþ becomes '-eþ' in Old English.")
    oe_form = "sēċeþ"
    print(f"Resulting Old English form: {oe_form}\n")
    
    # --- Step 4: Old High German ---
    print("--- 4. Derivation to Old High German (OHG) ---")
    print(f"Starting from PGmc: {pgmc_form}")
    print("  a) Diphthongization: The long vowel *ō becomes the diphthong 'uo'.")
    print("     *sōk- > suoh-")
    print("  b) High German Consonant Shift (part 1): The consonant *k becomes the fricative 'hh' (pronounced /x/) between vowels.")
    print("     The verb stem becomes 'suohh-' (as in the infinitive 'suohhen'). The 3sg form is 'suohhit'.")
    print("  c) High German Consonant Shift (part 2): The ending *-þ becomes '-t' for this verb class.")
    print("     *-iþ > -it")
    ohg_form = "suohhit"
    print(f"Resulting Old High German form: {ohg_form}\n")

    print("--- Summary ---")
    print(f"PIE:             *seh₂gieti")
    print(f"Proto-Germanic:  {pgmc_form}")
    print(f"Old Norse:       {on_form}")
    print(f"Old English:     {oe_form}")
    print(f"Old High German: {ohg_form}")


trace_etymology()