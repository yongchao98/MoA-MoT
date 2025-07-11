def derive_middle_english_word():
    """
    This function outlines the hypothetical linguistic derivation of a word
    from Proto-Indo-European to Middle English based on standard sound changes.
    """
    # Step 1: From Proto-Indo-European to Proto-Germanic
    step1_intro = "1. From Proto-Indo-European (PIE) to Proto-Germanic (PGmc):"
    pie_root = "*kʷeys-"
    pie_causative = "*kʷoyséye-"
    pgmc_form = "*hʷaizjaną"
    step1_desc = (
        f"   - Start with PIE root {pie_root} ('to see').\n"
        f"   - Form the o-grade causative to mean 'to show', resulting in {pie_causative}.\n"
        f"   - Apply Grimm's Law (*kʷ > *hʷ) and Verner's Law (*s > *z).\n"
        f"   - Adjust vowels (*oy > *ai) and suffix to get the PGmc infinitive: {pgmc_form}."
    )

    # Step 2: From Proto-Germanic to Old English
    step2_intro = "2. From Proto-Germanic (PGmc) to Old English (OE):"
    pgmc_infinitive = "*hʷaizjaną"
    oe_form = "hwǣrþ"
    step2_desc = (
        f"   - The /j/ in {pgmc_infinitive} triggers i-umlaut (*ai > ǣ) and rhotacism (*z > r).\n"
        f"   - This creates the OE infinitive 'hwǣran'.\n"
        f"   - The 3rd person singular present form adds the ending '-þ', giving: {oe_form}."
    )

    # Step 3: From Old English to Middle English
    step3_intro = "3. From Old English (OE) to Middle English (ME):"
    me_form = "whereth"
    step3_desc = (
        f"   - The OE vowel ǣ shifts to ME ē (spelled 'e').\n"
        f"   - The OE spelling 'hw' becomes 'wh'.\n"
        f"   - The OE ending '-þ' becomes the common ME ending '-eth'.\n"
        f"   - The final ME form for 'he shows' is: {me_form}."
    )

    # Print all steps and the final answer
    print("Derivation Steps:\n")
    print(step1_intro)
    print(step1_desc)
    print("\n" + step2_intro)
    print(step2_desc)
    print("\n" + step3_intro)
    print(step3_desc)
    print("\n-------------------------")
    print("Final Answer:")
    print(me_form)

derive_middle_english_word()
<<<whereth>>>