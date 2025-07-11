def derive_middle_english_verb():
    """
    Traces the hypothetical evolution of a PIE root to a Middle English verb form,
    printing each step of the derivation.
    """

    # --- Step 1: Proto-Indo-European (PIE) ---
    pie_root = "*kʷeys-"
    print(f"1. Start with the Proto-Indo-European root: {pie_root} (to see, to heed)")

    # --- Step 2: Form PIE O-Grade Causative ---
    pie_causative = "*kʷoys-eye-"
    print(f"2. Form the o-grade causative stem: {pie_causative} (to cause to see, i.e., 'to show')")

    # --- Step 3: Evolve to Proto-Germanic (PGmc) ---
    pgmc_infinitive = "*hʷazijaną"
    print(f"\n3. Apply sound changes to get the Proto-Germanic form:")
    print(f"   - Grimm's Law changes *kʷ > *hʷ")
    print(f"   - Verner's Law changes *s > *z (due to accent on the suffix)")
    print(f"   - Vowel shifts change *o > *a and the suffix *-eye- becomes *-ijaną")
    print(f"   - The resulting PGmc infinitive is: {pgmc_infinitive}")

    # --- Step 4: Evolve to Old English (OE) ---
    oe_form = "hwæreþ"
    print(f"\n4. Apply sound changes to get the Old English form:")
    print(f"   - I-Umlaut (the 'j' fronts the root vowel 'a' to 'æ')")
    print(f"   - Rhotacism (the 'z' becomes 'r')")
    print(f"   - Add the 3rd person singular present ending '-eþ'")
    print(f"   - The resulting OE 3rd p. sg. present form is: {oe_form}")

    # --- Step 5: Evolve to Middle English (ME) ---
    me_form = "whareth"
    print(f"\n5. Apply sound changes to get the Middle English form:")
    print(f"   - The spelling of 'hw' becomes 'wh'")
    print(f"   - The vowel 'æ' lowers to 'a'")
    print(f"   - The ending '-eþ' is spelled '-eth'")
    print(f"   - The final ME 3rd p. sg. present form is: {me_form}")

derive_middle_english_verb()