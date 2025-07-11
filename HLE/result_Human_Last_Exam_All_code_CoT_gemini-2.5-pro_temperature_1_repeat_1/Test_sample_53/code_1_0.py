def solve_linguistic_derivation():
    """
    This function traces the hypothetical evolution of a PIE root to a Middle English
    verbal form, printing each step of the process.
    """
    # Define the forms at each stage of the derivation
    pie_causative_3sg = "*kʷoyséyeti"
    pgmc_3sg = "*hʷaiziþi"
    oe_3sg = "hwǣreþ"
    me_3sg = "whereth"

    print("Derivation of a hypothetical Middle English verb from PIE *kʷeys- ('to see, to heed')")
    print("Target: 3rd person singular present, o-grade causative ('he shows')")
    print("-" * 75)

    # Step 1: From Proto-Indo-European to Proto-Germanic
    print("Step 1: Proto-Indo-European (PIE) -> Proto-Germanic (PGmc)")
    print(f"  - The PIE o-grade causative 3sg form is reconstructed as ~{pie_causative_3sg}.")
    print("  - Sound changes applied:")
    print("    - Grimm's Law: *kʷ > *hʷ")
    print("    - Verner's Law: *s > *z (due to non-initial-syllable accent in PIE)")
    print("    - Vowel Shift: *oy > *ai")
    print(f"  - Resulting PGmc stem is *hʷaiz-, giving the 3sg form ~{pgmc_3sg}.")
    print("-" * 75)

    # Step 2: From Proto-Germanic to Old English
    print("Step 2: Proto-Germanic (PGmc) -> Old English (OE)")
    print(f"  - The PGmc 3sg form is ~{pgmc_3sg}.")
    print("  - Sound changes applied:")
    print("    - i-Umlaut: The *i in the ending causes the root vowel *ai to become ǣ.")
    print("    - Rhotacism: *z becomes r.")
    print("    - The 3sg ending becomes -eþ.")
    print(f"  - Resulting OE form is '{oe_3sg}'.")
    print("-" * 75)

    # Step 3: From Old English to Middle English
    print("Step 3: Old English (OE) -> Middle English (ME)")
    print(f"  - The OE form is '{oe_sg}'.")
    print("  - Sound and spelling changes applied:")
    print("    - The initial 'hw' is spelled 'wh'.")
    print("    - The long vowel 'ǣ' becomes the long close vowel 'ē', spelled 'e'.")
    print("    - The ending '-eþ' is preserved as '-eth'.")
    print(f"  - Resulting ME form is '{me_3sg}'.")
    print("-" * 75)

    # Final Summary
    print("Final Derivation Equation:")
    # Using the stem for clarity in the early stages
    print("PIE *kʷoys- > PGmc *hʷaiz- > OE hwǣr- > ME wher-")
    print("\nFull 3rd Person Singular Present Verb Form Evolution:")
    print(f"{pie_causative_3sg} > {pgmc_3sg} > {oe_3sg} > {me_3sg}")
    print("-" * 75)

    print("\nThe final predicted Middle English form is:")
    print(me_3sg)


solve_linguistic_derivation()
<<<whereth>>>