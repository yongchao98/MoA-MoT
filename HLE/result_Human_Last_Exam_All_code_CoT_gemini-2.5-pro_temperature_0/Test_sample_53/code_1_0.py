def trace_word_evolution():
    """
    Traces a hypothetical PIE root through its sound changes to Middle English.
    """
    # Step 1: Proto-Indo-European (PIE)
    # The task is to find the reflex of the PIE root *kʷeys- ('to see, to heed')
    # as an o-grade causative verb in the 3rd person singular present.
    pie_root = "*kʷeys-"
    # The o-grade changes the root vowel 'e' to 'o' -> *kʷoys-
    # The causative is formed with the suffix -éye-.
    # The 3rd person singular present ending is -ti.
    pie_form = "*kʷoyséyeti"
    print("1. Proto-Indo-European Stage:")
    print(f"   - PIE Root: {pie_root}")
    print(f"   - O-grade causative 3sg. present form: {pie_form}")
    print("-" * 20)

    # Step 2: Proto-Germanic (PGmc)
    # We apply sound changes from PIE to PGmc.
    # - Grimm's Law: *kʷ -> *hʷ and *t -> *þ
    # - Verner's Law: *s -> *z (because the PIE accent was on the suffix, not the root)
    # - Vowel/Suffix changes: PIE *oy -> PGmc *ai; the causative suffix *-éye- becomes *-i- in this form.
    pgmc_form = "*hwaiziþi"
    print("2. Proto-Germanic Stage:")
    print(f"   - Sound changes applied: Grimm's Law, Verner's Law, vowel shifts.")
    print(f"   - The form becomes: {pgmc_form}")
    print("-" * 20)

    # Step 3: Old English (OE)
    # We apply sound changes from PGmc to OE.
    # - i-umlaut: The 'i' in the stem fronts the root diphthong *ai -> ǣ.
    # - Rhotacism: The voiced fricative *z becomes 'r'.
    # - Inflectional ending: The 3sg. ending *-iþ(i) becomes -eþ.
    # - Spelling: *hƿ is written as 'hƿ'.
    oe_form = "hƿǣreþ"
    print("3. Old English Stage:")
    print(f"   - Sound changes applied: i-umlaut, rhotacism, inflectional changes.")
    print(f"   - The form becomes: {oe_form}")
    print("-" * 20)

    # Step 4: Middle English (ME)
    # We apply sound changes from OE to ME.
    # - Spelling: 'hƿ' becomes 'wh'.
    # - Vowel Shift: OE long 'ǣ' /æː/ becomes ME long open 'ē' /ɛː/, spelled 'e'.
    # - Inflectional ending: The ending -eþ is preserved as -eth.
    me_form = "whereth"
    print("4. Middle English Stage:")
    print(f"   - Sound changes applied: Spelling modernization, vowel shifts.")
    print(f"   - The final 3rd person singular present form is: {me_form}")
    print("-" * 20)
    
    print(f"The final predicted Middle English form for 'he shows' is '{me_form}'.")

trace_word_evolution()
<<<whereth>>>