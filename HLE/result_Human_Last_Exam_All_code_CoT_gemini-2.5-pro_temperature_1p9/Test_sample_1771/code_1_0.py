def trace_germanic_etymology():
    """
    Traces the etymological development of a PIE verb into Germanic languages.
    The script simulates major sound changes from a reconstructed Proto-Germanic form.
    """

    # Step 1: Establish the Proto-Germanic form.
    # The PIE verb *seh₂gieti ('(s)he gives a sign, tracks') gives rise to the
    # Proto-Germanic (PGmc) weak verb *sōkijaną ('to seek'). We use the
    # 3rd person singular present indicative form, reconstructed as *sōkiþi,
    # to match the PIE ending *-ti.
    pgmc_form = "*sōkiþi"
    pgmc_root = "sōk"
    pgmc_ending_sound = "iþi"
    
    print(f"The PIE verb *seh₂gieti evolves into the following forms in the Germanic languages.")
    print("-" * 20)
    print(f"Proto-Germanic: {pgmc_form} (root: '{pgmc_root}-' + ending: '-{pgmc_ending_sound}')")
    print("-" * 20)

    # Step 2: Derive the Old Norse form
    # Rule 1 (i-umlaut): The back vowel 'ō' is fronted to 'œ' by the following 'i'.
    on_root = pgmc_root.replace("ō", "œ")
    # Rule 2 (Ending): The PGmc ending *-iþi evolves into Old Norse '-ir' for this verb class.
    on_ending = "ir"
    on_form = on_root + on_ending
    print("Old Norse Derivation:")
    print(f"1. Root '{pgmc_root}-' undergoes i-umlaut because of the following 'i' -> '{on_root}-'")
    print(f"2. Ending '-{pgmc_ending_sound}' becomes -> '-{on_ending}'")
    print(f"Final Form: '{on_root}-' + '-{on_ending}' = {on_form}")
    print("-" * 20)

    # Step 3: Derive the Old English form
    # Rule 1 (i-umlaut/Anglo-Frisian Brightening): 'ō' becomes 'ē'.
    oe_root_step1 = pgmc_root.replace("ō", "ē")
    # Rule 2 (Palatalization): 'k' becomes 'ċ' (pronounced /tʃ/) before a front vowel like 'ē' or 'i'.
    oe_root_step2 = oe_root_step1.replace("k", "ċ")
    # Rule 3 (Ending): The PGmc ending *-iþi evolves into Old English '-eþ'.
    oe_ending = "eþ"
    oe_form = oe_root_step2 + oe_ending
    print("Old English Derivation:")
    print(f"1. Root '{pgmc_root}-' undergoes i-umlaut -> '{oe_root_step1}-'")
    print(f"2. Root '{oe_root_step1}-' undergoes palatalization -> '{oe_root_step2}-'")
    print(f"3. Ending '-{pgmc_ending_sound}' becomes -> '-{oe_ending}'")
    print(f"Final Form: '{oe_root_step2}-' + '-{oe_ending}' = {oe_form}")
    print("-" * 20)

    # Step 4: Derive the Old High German form
    # Rule 1 (Diphthongization): 'ō' becomes 'uo'.
    ohg_root_step1 = pgmc_root.replace("ō", "uo")
    # Rule 2 (High German Consonant Shift on root): 'k' after a vowel becomes 'hh' (a long fricative /xx/).
    ohg_root_step2 = ohg_root_step1.replace("k", "hh")
    # Rule 3 (High German Consonant Shift on ending): 'þ' becomes 'd'. The ending *-iþi -> *-idi -> -it.
    ohg_ending = "it"
    ohg_form = ohg_root_step2 + ohg_ending
    print("Old High German Derivation:")
    print(f"1. Root '{pgmc_root}-' undergoes diphthongization -> '{ohg_root_step1}-'")
    print(f"2. Root '{ohg_root_step1}-' undergoes consonant shift -> '{ohg_root_step2}-'")
    print(f"3. Ending '-{pgmc_ending_sound}' becomes -> '-{ohg_ending}'")
    print(f"Final Form: '{ohg_root_step2}-' + '-{ohg_ending}' = {ohg_form}")
    print("-" * 20)

# Run the tracing function
trace_germanic_etymology()
<<<Proto-Germanic: *sōkiþi, Old Norse: sœkir, Old English: sēċeþ, Old High German: suohhit>>>