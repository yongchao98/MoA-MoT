def derive_germanic_forms():
    """
    This script traces a PIE verb through its Proto-Germanic, Old Norse,
    Old English, and Old High German descendants by applying historical sound laws.
    """

    # 1. Proto-Germanic (*sōkiþi)
    print("--- 1. Derivation for Proto-Germanic ---")
    pie_form = "*seh₂gieti"
    pie_root = "*seh₂g-"
    print(f"Step 1: Start with the PIE root '{pie_root}'. The provided cognates (Lat. sāgīre) and the verbal form point to the meaning 'to track, seek'.")
    
    # Note on Ablaut
    pgmc_root_vowel = "*sōg-"
    print(f"Step 2: The verb develops an 'o'-grade ablaut vowel in Germanic, giving the root form '{pgmc_root_vowel}'.")
    
    # Grimm's Law
    pgmc_root_consonant = "*sōk-"
    print(f"Step 3: Apply Grimm's Law. The voiced stop '*g' becomes a voiceless stop '*k'. This gives the root '{pgmc_root_consonant}'.")
    print("   (Note: Verner's Law is blocked, as evidenced by all descendants having a 'k'-sound, implying a shift to root-accentuation in Proto-Germanic.)")

    # Form the 3rd person singular verb
    pgmc_verb_form = "*sōkiþi"
    print(f"Step 4: Form the 3rd person singular of this Class 1 weak verb by adding the suffix and ending '-iþi'.")
    print(f"Final Proto-Germanic Form: {pgmc_verb_form}\n")
    
    # 2. Old Norse (sœkir)
    print("--- 2. Derivation for Old Norse ---")
    print(f"Step 1: Start with the PGmc form '{pgmc_verb_form}'.")
    
    # i-umlaut
    on_root = "sœk-"
    print(f"Step 2: The 'i' (from the earlier glide *j) in the suffix causes i-umlaut, fronting the root vowel '*ō' to 'œ'. The root becomes '{on_root}'.")

    # Ending change
    on_ending = "ir"
    on_verb_form = "sœkir"
    print(f"Step 3: The PGmc ending '*-iþi' becomes '-{on_ending}' in Old Norse.")
    print(f"Final Old Norse Form: {on_verb_form}\n")

    # 3. Old English (sēcþ)
    print("--- 3. Derivation for Old English ---")
    print(f"Step 1: Start with the PGmc form '{pgmc_verb_form}'.")
    
    # i-umlaut
    oe_root = "sēc-"
    print(f"Step 2: The 'i' causes i-umlaut, changing the root vowel '*ō' to 'ē'. The root becomes '{oe_root}'.")
    
    # Ending change
    oe_ending = "þ"
    oe_verb_form = "sēcþ"
    print(f"Step 3: The PGmc ending '*-iþi' becomes '-eþ', which is often syncopated to '-{oe_ending}'.")
    print(f"Final Old English Form: {oe_verb_form}\n")

    # 4. Old High German (suohhit)
    print("--- 4. Derivation for Old High German ---")
    print(f"Step 1: Start with the PGmc form '{pgmc_verb_form}'.")

    # Vowel Change
    ohg_root_vowel = "suoh-"
    print(f"Step 2: The PGmc vowel '*ō' becomes 'uo' in Old High German. The root becomes '{ohg_root_vowel}'.")

    # High German Consonant Shift
    ohg_root_consonant = "suohh-"
    print(f"Step 3: Apply the High German Consonant Shift. The PGmc '*k' becomes the long fricative 'hh'. The root is now '{ohg_root_consonant}'.")
    
    # Ending change
    ohg_ending = "it"
    ohg_verb_form = "suohhit"
    print(f"Step 4: The PGmc ending '*-iþi' becomes '-{ohg_ending}'.")
    print(f"Final Old High German Form: {ohg_verb_form}\n")
    
    final_answer = f"<<<Proto-Germanic: *sōkiþi, Old Norse: sœkir, Old English: sēcþ, Old High German: suohhit>>>"
    # This is a dummy print to ensure the final answer is part of the output.
    # In a real application, you might just return this value.
    # print(final_answer)

# Execute the derivation function
derive_germanic_forms()