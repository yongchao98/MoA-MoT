def trace_etymology():
    """
    Traces the PIE form *seh₂gieti through Germanic languages
    and prints the derivation.
    """
    # 1. Define the forms at each stage
    pie_form = "*seh₂gieti"
    
    # Proto-Germanic
    # Grimm's Law: *g > *k
    # Vowel change: *eh₂ > *ā > *ō
    # Ending change: *-ieti > *-iþ (3sg. pres. weak verb)
    pgmc_form = "*sōkiþ"
    
    # Old Norse from PGmc *sōkiþ
    # i-umlaut: *ō > œ
    # Ending: *-iþ > -ir
    on_form = "sœkir"
    
    # Old English from PGmc *sōkiþ
    # i-umlaut: *ō > ē
    # Palatalization: *k > ċ
    # Ending: *-iþ > -eþ
    oe_form = "sēċeþ"
    
    # Old High German from PGmc *sōkiþ
    # Diphthongization: *ō > uo
    # High German Consonant Shift: *k > hh
    # Ending: *-iþ > -it
    ohg_form = "suohhit"

    # 2. Print the explanation and results
    print(f"The derivation of PIE '{pie_form}' is as follows:\n")

    print("1. Proto-Germanic:")
    print(f"   The PIE form becomes '{pgmc_form}'.")
    print("   - Grimm's Law changes the root consonant *g -> *k.")
    print("   - The PIE vowel and laryngeal *eh₂ become the long vowel *ō.")
    print("   - The 3rd person singular ending becomes *-iþ.\n")

    print("2. Old Norse:")
    print(f"   Proto-Germanic '{pgmc_form}' becomes '{on_form}'.")
    print("   - The root vowel *ō undergoes i-umlaut to become œ.")
    print("   - The ending *-þ becomes -r.\n")

    print("3. Old English:")
    print(f"   Proto-Germanic '{pgmc_form}' becomes '{oe_form}'.")
    print("   - The root vowel *ō undergoes i-umlaut to become ē.")
    print("   - The consonant *k is palatalized to ċ (/tʃ/) by the following 'i'.")
    print("   - The ending *-iþ becomes -eþ.\n")

    print("4. Old High German:")
    print(f"   Proto-Germanic '{pgmc_form}' becomes '{ohg_form}'.")
    print("   - The root vowel *ō becomes the diphthong uo.")
    print("   - The High German Consonant Shift changes *k to hh (/x/).")
    print("   - The ending *-iþ becomes -it.\n")
    
    print("--- Final Derivation Chain ---")
    # The final "equation" with each form outputted
    print(f"{pie_form} > {pgmc_form} > ON {on_form}, OE {oe_form}, OHG {ohg_form}")

trace_etymology()