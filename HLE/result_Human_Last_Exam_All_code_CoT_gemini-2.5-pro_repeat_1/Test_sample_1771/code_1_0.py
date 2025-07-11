def trace_etymology():
    """
    Traces a PIE word through Proto-Germanic, Old Norse, Old English,
    and Old High German, explaining the sound changes at each step.
    """
    pie_form = "*seh₂giéti"
    
    print(f"Tracing the development of PIE '{pie_form}' ('(s)he seeks/tracks').\n")

    # Step 1: Proto-Germanic
    print("--- 1. Proto-Germanic ---")
    print("a. The root is *seh₂g-. The laryngeal *h₂ causes the preceding *e to become *a, giving *sāg-.")
    print("b. In Proto-Germanic, PIE long *ā regularly becomes *ō. The root is now *sōg-.")
    print("c. By Grimm's Law, the PIE voiced stop *g becomes a voiceless stop *k. The root is now *sōk-.")
    print("d. The PIE 3rd person singular present ending *-ieti becomes *-iþi in Proto-Germanic for this verb class.")
    pgmc_form = "*sōkiþi"
    print(f"Result for Proto-Germanic: {pgmc_form}\n")

    # Step 2: Old Norse
    print("--- 2. Old Norse (from PGmc *sōkiþi) ---")
    print("a. The *i in the syllable after the root causes 'i-umlaut', fronting the root vowel *ō to œ.")
    print("b. The ending *-iþi evolves into the 3rd person singular ending -ir.")
    on_form = "sœkir"
    print(f"Result for Old Norse: {on_form}\n")

    # Step 3: Old English
    print("--- 3. Old English (from PGmc *sōkiþi) ---")
    print("a. 'i-umlaut' also occurs here, fronting the root vowel *ō to ē.")
    print("b. The consonant *k is palatalized to /tʃ/ (written 'c') before a front vowel like the 'i' that followed.")
    print("c. The ending *-iþi becomes -eþ in Old English.")
    oe_form = "sēceþ"
    print(f"Result for Old English: {oe_form}\n")

    # Step 4: Old High German
    print("--- 4. Old High German (from PGmc *sōkiþi) ---")
    print("a. The root vowel *ō 'breaks' into the diphthong 'uo'.")
    print("b. The High German Consonant Shift changes *k to the fricative 'hh' (often written 'h') and *þ to 'd'.")
    print("c. The ending evolves into -it. Combining these gives *suohh-id-it -> suohhit.")
    ohg_form = "suohit"
    print(f"Result for Old High German: {ohg_form}\n")

    # Final Summary
    final_answer = (f"Proto-Germanic: {pgmc_form}, Old Norse: {on_form}, "
                    f"Old English: {oe_form}, Old High German: {ohg_form}")
    
    # The final answer is wrapped according to the format requirements.
    print("--- FINAL SUMMARY ---")
    print(final_answer)
    print(f"\n<<<{final_answer}>>>")

# Execute the function
trace_etymology()