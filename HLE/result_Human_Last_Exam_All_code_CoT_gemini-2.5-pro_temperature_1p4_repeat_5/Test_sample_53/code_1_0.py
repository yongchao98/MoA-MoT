def derive_hypothetical_word():
    """
    Traces a hypothetical PIE root through its sound changes into Middle English.
    """
    # Step 1: Define the Proto-Indo-European (PIE) root.
    pie_root = "*kʷeys-"
    print(f"1. Proto-Indo-European Root: {pie_root} ('to see, to heed')")

    # Step 2: Form the o-grade causative in PIE. This means "to cause to see/heed".
    pie_causative_stem = "*kʷoyséye-"
    print(f"2. PIE O-Grade Causative Stem: {pie_causative_stem}")

    # Step 3: Apply sound changes to get the Proto-Germanic (PGmc) form.
    # Grimm's Law (*kʷ > *hʷ), Verner's Law (*s > *z), Vowel shift (*oy > *ai), Suffix change (*-eye- > *-ijaną)
    pgmc_infinitive = "*hʷaizijaną"
    print(f"3. Proto-Germanic Infinitive: {pgmc_infinitive}")

    # Step 4: Trace to Old English (OE) 3rd person singular present.
    # i-Umlaut (*ai > ǣ), Rhotacism (*z > r), and addition of the '-eþ' ending.
    old_english_form = "hwǣreþ"
    print(f"4. Old English 3rd Person Singular Present: {old_english_form}")

    # Step 5: Trace to Middle English (ME) 3rd person singular present.
    # Vowel shift (ǣ > ē, spelled 'e'), spelling changes (hw > wh, þ > th).
    middle_english_form = "whereth"
    print(f"5. Middle English 3rd Person Singular Present: {middle_english_form} ('he shows')")

# Execute the function to show the derivation.
derive_hypothetical_word()