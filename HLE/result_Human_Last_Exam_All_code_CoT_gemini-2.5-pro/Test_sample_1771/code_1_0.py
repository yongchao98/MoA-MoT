def trace_linguistic_evolution():
    """
    Traces the evolution of a PIE word through Germanic languages
    by printing the application of sound laws.
    """
    pie = "*seh₂gieti"
    print(f"Starting with the Proto-Indo-European (PIE) form: {pie}\n")

    # Step 1: Proto-Germanic
    print("--- 1. Deriving Proto-Germanic (PGmc) ---")
    pgmc_step1 = "*sāg-"
    print(f"Laryngeal coloring: *eh₂ becomes *ā. The root becomes {pgmc_step1}.")
    pgmc_step2 = "*sōg-"
    print(f"Vowel shift: PIE *ā becomes PGmc *ō. The root becomes {pgmc_step2}.")
    pgmc_step3 = "*sōk-"
    print(f"Grimm's Law: The voiced stop *g becomes the voiceless stop *k. The root is now {pgmc_step3}.")
    pgmc_final = "*sōkijaną"
    print(f"The verb belongs to the Class 1 weak verbs, giving the infinitive form: {pgmc_final}.\n")

    # Step 2: Old Norse
    print("--- 2. Deriving Old Norse (ON) from PGmc *sōkijaną ---")
    on_result = "sœkja"
    print(f"I-umlaut: The *j fronts the root vowel *ō to œ.")
    print(f"Result: The Old Norse form is {on_result}.\n")

    # Step 3: Old English
    print("--- 3. Deriving Old English (OE) from PGmc *sōkijaną ---")
    oe_result = "sēcan"
    print(f"I-umlaut: The *j fronts the root vowel *ō to ē.")
    print(f"Palatalization: *k becomes ċ ([tʃ]) before the new front vowel ē.")
    print(f"Result: The Old English form is {oe_result}.\n")

    # Step 4: Old High German
    print("--- 4. Deriving Old High German (OHG) from PGmc *sōkijaną ---")
    ohg_result = "suohhen"
    print(f"High German Consonant Shift: *k becomes hh ([x]) after a vowel.")
    print(f"Vowel Shift: The root vowel *ō becomes the diphthong uo.")
    print(f"Result: The Old High German form is {ohg_result}.\n")

    # Final summary
    print("--- Summary of results ---")
    print(f"Proto-Germanic: {pgmc_final}")
    print(f"Old Norse: {on_result}")
    print(f"Old English: {oe_result}")
    print(f"Old High German: {ohg_result}")

if __name__ == "__main__":
    trace_linguistic_evolution()