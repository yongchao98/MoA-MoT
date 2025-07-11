def derive_word():
    """
    This function traces the hypothetical derivation of a word from Proto-Indo-European
    to Middle English and prints the result at each stage, explaining the sound laws applied.
    """
    
    # Stage 1: Proto-Indo-European (PIE)
    # The root is *kʷeys- (to see, to heed).
    # We form the o-grade causative, 3rd person singular present.
    # The o-grade of *kʷeys- is *kʷoys-.
    # The 3rd person singular causative form is *kʷoyséyeti.
    pie_form = "*kʷoyséyeti"
    print(f"1. Proto-Indo-European form: {pie_form}")
    print("   - Explanation: The root *kʷeys- is put into the o-grade (*kʷoys-) and combined with the causative suffix and 3rd person singular ending.")
    print("-" * 30)

    # Stage 2: Proto-Germanic (PGmc)
    # Changes: Grimm's Law (*kʷ > *hʷ), Verner's Law (*s > *z), vowel shift (*oy > *ai), suffix simplification.
    pgmc_form = "*hwaiziþ"
    print(f"2. Proto-Germanic form: {pgmc_form}")
    print("   - Explanation: Grimm's Law changes *kʷ to *hʷ. Verner's Law changes *s to *z due to the unaccented root syllable in PIE. The diphthong *oy becomes *ai, and the verbal ending simplifies to *-iþ for a Class 1 weak verb.")
    print("-" * 30)

    # Stage 3: Old English (OE)
    # Changes: i-umlaut (*ai > ǣ), rhotacism (*z > r), syncope.
    # The infinitive PGmc *hwaizjaną shows the conditions for umlaut and rhotacism, which then apply to the whole verb pattern.
    oe_form = "hwǣrþ"
    print(f"3. Old English form: {oe_form}")
    print("   - Explanation: The stem vowel *ai undergoes i-umlaut to ǣ. The consonant *z undergoes rhotacism to become r. The vowel in the ending *-iþ is lost (syncope) after a long root syllable, leaving -þ.")
    print("-" * 30)
    
    # Stage 4: Middle English (ME)
    # Changes: spelling change (hw- > wh-), vowel shift (ǣ > ē), ending change (-þ > -eth).
    me_form = "whereth"
    print(f"4. Middle English form: {me_form}")
    print("   - Explanation: The OE cluster hw- is now spelled wh-. The long vowel ǣ raises to a long open ē (often spelled 'e' or 'ee'). The 3rd person singular ending -þ becomes -eth in many dialects.")
    print("-" * 30)

    print(f"\nThe final derived form in Middle English is '{me_form}'.")


derive_word()
print("\n<<<whereth>>>")