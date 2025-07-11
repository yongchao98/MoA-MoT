def derive_middle_english_verb():
    """
    This function traces the hypothetical phonological history of a PIE root
    to a Middle English verb form, explaining each step of the derivation.
    """
    
    # Define the initial forms and their meanings.
    pie_root = "*kʷeys-"
    pie_meaning = "to see, to heed"
    derivation_type = "o-grade causative"
    
    print("This script will derive a hypothetical Middle English word from a Proto-Indo-European root.")
    print("-" * 50)

    # Step 1: From Proto-Indo-European to Proto-Germanic
    print("Step 1: Proto-Indo-European -> Proto-Germanic\n")
    pie_causative_3sg = "*kʷoyséyeti"
    pgmc_infinitive = "*hwaizijaną"
    print(f"We begin with the PIE root {pie_root} ({pie_meaning}).")
    print(f"The o-grade causative 3rd person singular form would be approximately {pie_causative_3sg}.")
    print("\nApplying sound changes to get to Proto-Germanic:")
    print("1. Grimm's Law: *kʷ > *hw")
    print("2. Verner's Law: *s > *z (due to unstressed root syllable in Proto-Germanic)")
    print("3. Vowel Shift: PIE *o > PGmc *a")
    print("4. Suffix Change: The causative becomes a Class 1 Weak Verb, with infinitive suffix *-ijaną.")
    print(f"\nThe resulting Proto-Germanic infinitive is: {pgmc_infinitive} (meaning 'to show' or 'to make aware').")
    print("-" * 50)
    
    # Step 2: From Proto-Germanic to Old English
    print("Step 2: Proto-Germanic -> Old English\n")
    pwgmc_form = "*hwairijan"
    oe_3sg_present = "hwǣreþ"
    print("First, in the transition to Proto-West Germanic, rhotacism changes *z to *r:")
    print(f"  {pgmc_infinitive} > {pwgmc_form}\n")
    print("Next, applying changes to get to Old English:")
    print("1. i-Umlaut: The *j in the suffix fronts the root diphthong *ai to ǣ.")
    print("2. The medial *j is then lost after the long vowel.")
    print(f"This results in the Old English infinitive 'hwǣran'.")
    print("3. For the 3rd person singular present, the ending '-eþ' is added to the stem 'hwǣr-'.")
    print(f"\nThe resulting Old English form is: {oe_3sg_present}")
    print("-" * 50)

    # Step 3: From Old English to Middle English
    print("Step 3: Old English -> Middle English\n")
    me_form = "whereth"
    print("Applying sound changes to get to Middle English:")
    print("1. Spelling: OE 'hw-' is regularly spelled 'wh-' in Middle English.")
    print("2. Vowel Shift: OE long 'ǣ' becomes the long open vowel 'ē' (written 'e').")
    print("3. Inflection: The verb ending OE '-eþ' becomes '-eth' in Midland/Southern Middle English.")
    print(f"\nThe resulting Middle English 3rd person singular present form is: {me_form}")
    print("-" * 50)

    # Final Summary
    print("Summary of the derivation chain:")
    # We print each of the reconstructed forms in the "equation" as requested
    chain_elements = [pie_causative_3sg, pgmc_infinitive, oe_3sg_present, me_form]
    print(f"{' > '.join(chain_elements)}")
    print("\nThe word means approximately 'he shows' or 'he causes (one) to heed'.")

derive_middle_english_verb()
<<<whereth>>>