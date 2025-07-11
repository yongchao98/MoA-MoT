def derive_middle_english_verb():
    """
    Prints the step-by-step derivation of a hypothetical Middle English verb
    from a Proto-Indo-European root, following standard sound changes.
    The instruction 'output each number in the final equation' is interpreted
    as printing each numbered step in the derivation.
    """
    derivation_steps = [
        ("Proto-Indo-European Root", "*kʷeys- (to see, to heed)"),
        ("1. PIE O-Grade Causative Formation", 
         "The o-grade of *kʷeys- is *kʷoys-. The 3rd person singular causative is formed with the suffix *-éye-ti, yielding *kʷoyséyeti, meaning 'he causes to see' or 'he shows'."),
        ("2. Transformation to Proto-Germanic (PGmc)", 
         "Applying sound laws:\n"
         "  - Grimm's Law: *kʷ > *hʷ.\n"
         "  - Verner's Law: Since the PIE accent is on the suffix (*-éye-), the unaccented root is followed by *s, which voices to *z.\n"
         "  - Vowel Shift: PIE *o > PGmc *a.\n"
         "  - Result: The PGmc verb stem is *hʷazjan- (a Class 1 weak verb). The 3rd person singular is *hʷaziþi."),
        ("3. Transformation to Old English (OE)", 
         "Applying sound laws from PGmc to OE:\n"
         "  - I-Umlaut: The 'i' in the PGmc form *hʷaziþi causes the root vowel 'a' to front to 'æ'.\n"
         "  - Consonant Development: PGmc *z becomes 'r' in Old English. PGmc *hʷ becomes 'hw'.\n"
         "  - Inflection: The 3rd person singular ending becomes '-eþ'.\n"
         "  - Result: The Old English form is hwæreþ."),
        ("4. Transformation to Middle English (ME)", 
         "Applying sound laws from OE to ME:\n"
         "  - Vowel Shift: Old English 'æ' becomes 'e'.\n"
         "  - Spelling Change: 'hw' is commonly spelled 'wh'.\n"
         "  - Inflection: The ending '-eþ' is spelled '-eth'.\n"
         "  - Final Result: The stem hwær- + eth becomes hwereth.")
    ]

    print("Deriving the Middle English form of PIE *kʷeys- ('he shows'):")
    print("-" * 60)
    for step, description in derivation_steps:
        print(f"Step {step}:\n{description}\n")

    # The "equation" is the derivation itself. Each part is printed above.
    # The final word represents the solution to this derivation.
    final_form = "hwereth"
    print(f"The final equation of sound changes results in the word: {final_form}")


if __name__ == '__main__':
    derive_middle_english_verb()
    print("\n<<<hwereth>>>")