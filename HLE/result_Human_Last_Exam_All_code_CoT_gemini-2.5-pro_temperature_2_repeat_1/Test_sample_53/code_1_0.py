import textwrap

def print_derivation():
    """
    Prints the step-by-step linguistic derivation of a hypothetical
    Middle English word from a Proto-Indo-European root.
    """
    # Introduction
    print("Deriving the Middle English reflex for PIE *kʷeys- ('to see, to heed').")
    print("Path: o-grade causative -> Proto-Germanic -> Old English -> Middle English 3rd Person Singular.\n")

    # Step 1: PIE to Proto-Germanic
    print("Step 1: Proto-Indo-European to Proto-Germanic")
    print("-" * 50)
    print(textwrap.fill("Initial PIE root: *kʷeys-", 70))
    print(textwrap.fill("Formation of o-grade causative stem: *kʷoys-eye- (meaning 'to cause to see', i.e., 'to show')", 70))
    print("\nSound Changes into Proto-Germanic:")
    print("  - Grimm's Law: *kʷ > *hʷ")
    print("  - Verner's Law (on *s, due to unaccented root): *s > *z")
    print("  - Vowel Shift: PIE *oy > PGmc *ai")
    print(f"Resulting PGmc verb stem: *hʷaiz-")
    print("Formation as a Class 1 Weak Verb (causative):")
    final_pgmc = "*hʷaizijaną"
    print(f"  Proto-Germanic Form: *hʷaizijaną (he shows = *hʷaiziþi)\n")


    # Step 2: Proto-Germanic to Old English
    print("Step 2: Proto-Germanic to Old English")
    print("-" * 50)
    print(f"Starting from PGmc infinitive: {final_pgmc}")
    print("\nSound Changes into Old English:")
    print("  - I-Umlaut: The 'j' in *-ijaną fronts the root vowel *ai > ǣ.")
    print("  - Rhotacism: Intervocalic *z > r.")
    print("  - Suffix evolution to standard infinitive '-an'.")
    print("Resulting Old English infinitive: hƿǣran (analogy: PGmc *laizijaną > OE lǣran)")
    final_oe_3ps = "hƿǣreþ"
    print(f"Old English 3rd person singular present form (stem hƿǣr- + suffix -eþ): {final_oe_3ps}\n")


    # Step 3: Old English to Middle English
    print("Step 3: Old English to Middle English")
    print("-" * 50)
    print(f"Starting from Old English 3rd person singular: {final_oe_3ps}")
    print("\nSound Changes into Middle English:")
    print("  - Spelling/Consonant Shift: OE 'hƿ' > ME 'wh'.")
    print("  - Vowel Shift: OE 'ǣ' (from i-mutated 'ai') > ME 'e' [ɛː].")
    print("  - Inflectional Ending Shift: OE '-eþ' > ME '-eth'.")
    print("\nCombining the changes:")
    print("  - Stem: hƿǣr- > wher-")
    print("  - Ending: -eþ > -eth")
    final_me_word = "whereth"
    print("\nFinal reconstructed Middle English form (meaning 'he shows'):")
    print(f"wh + er + eth = {final_me_word}")

print_derivation()