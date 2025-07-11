def trace_linguistic_evolution():
    """
    This function traces the PIE verb *seh2gieti to its descendants
    in several old Germanic languages and prints the derivation.
    """
    
    # The initial PIE form provided by the user, including the number 2.
    pie_form_original = "*seh2gieti"
    
    # We start from the o-grade of the root, *soh₂g-, which is the basis for the Germanic verb.
    pie_root_ograde = "*soh₂g-"
    
    # The derived forms in each language stage.
    proto_germanic = "*sōkijaną"
    old_norse = "sœkja"
    old_english = "sēċan"
    old_high_german = "suohhen"
    
    print(f"Tracing the evolution of PIE {pie_form_original} (root: {pie_root_ograde}):\n")
    
    # --- Step 1: Proto-Germanic ---
    print("1. To Proto-Germanic:")
    print(f"   PIE '{pie_root_ograde}' becomes PGmc. '{proto_germanic}'")
    print("   - Grimm's Law: The voiced stop *g shifts to the voiceless stop *k.")
    print("   - Vowel Change: The o-grade root with a laryngeal (*soh₂g-) yields a long vowel *ō.")
    print("   - Suffix: The PIE suffix evolves into the PGmc. infinitive *-janą.")
    print(f"   Result: {proto_germanic} (to seek)\n")

    # --- Step 2: Old Norse ---
    print("2. To Old Norse:")
    print(f"   PGmc. '{proto_germanic}' becomes ON '{old_norse}'")
    print("   - I-Umlaut: The *j in the suffix causes the root vowel *ō to be fronted to œ.")
    print("   - Ending: The infinitive ending *-aną reduces to -a.")
    print(f"   Result: {old_norse}\n")

    # --- Step 3: Old English ---
    print("3. To Old English:")
    print(f"   PGmc. '{proto_germanic}' becomes OE '{old_english}'")
    print("   - I-Umlaut: The *j causes the root vowel *ō to be fronted to ē.")
    print("   - Palatalization: The consonant *k is palatalized to ċ (pronounced /tʃ/) before the *j.")
    print("   - Suffix: The *j is lost and the infinitive ending becomes -an.")
    print(f"   Result: {old_english}\n")

    # --- Step 4: Old High German ---
    print("4. To Old High German:")
    print(f"   PGmc. '{proto_germanic}' becomes OHG '{old_high_german}'")
    print("   - High German Consonant Shift: The voiceless stop *k shifts to the fricative hh (written ch).")
    print("   - Vowel Change: The long vowel *ō diphthongizes to uo.")
    print("   - Suffix: The infinitive ending becomes -en.")
    print(f"   Result: {old_high_german}\n")

# Execute the function to print the results.
trace_linguistic_evolution()