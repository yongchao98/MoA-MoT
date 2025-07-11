def trace_word_derivation():
    """
    Traces the hypothetical derivation of a PIE root to a Middle English verb form.
    This script explains the sound changes at each stage of the derivation.
    """

    # --- Initial State: Proto-Indo-European (PIE) ---
    pie_root = "*kʷeys-"
    pie_causative = "*kʷoyséye-"
    print(f"1. Proto-Indo-European Stage")
    print(f"   - Start with the PIE root: {pie_root} ('to see, to heed')")
    print(f"   - Form the o-grade causative verb stem: {pie_causative} ('he/she causes to see')")
    print("-" * 20)

    # --- Step 2: PIE to Proto-Germanic (PGmc) ---
    pgmc_stem = "*hʷaizi-"
    pgmc_infinitive = "*hʷaizijaną"
    print(f"2. Proto-Germanic Stage (via sound changes)")
    print(f"   - Grimm's Law: PIE *kʷ -> PGmc *hʷ")
    print(f"   - Verner's Law: PIE *s -> PGmc *z (due to unaccented root)")
    print(f"   - Vowel Shift: PIE *oy -> PGmc *ai")
    print(f"   - The resulting PGmc verb is a Class 1 Weak Verb: {pgmc_infinitive}")
    print(f"   - The effective stem is: {pgmc_stem}")
    print("-" * 20)

    # --- Step 3: PGmc to Old English (OE) ---
    oe_infinitive = "hwǣran"
    oe_3sg_present = "hwǣreþ"
    print(f"3. Old English Stage")
    print(f"   - I-mutation: The '-ij-' in the PGmc infinitive suffix mutates the root vowel *ai -> ǣ")
    print(f"   - Rhotacism: Medial *z between vowels becomes 'r' in Old English")
    print(f"   - Initial Consonant: PGmc *hʷ becomes OE 'hw'")
    print(f"   - The resulting OE infinitive is: {oe_infinitive}")
    print(f"   - The 3rd person singular present form adds the ending '-eþ': {oe_3sg_present}")
    print("-" * 20)

    # --- Step 4: OE to Middle English (ME) ---
    me_infinitive = "wheren"
    me_3sg_present = "whereth"
    print(f"4. Middle English Stage")
    print(f"   - Vowel Shift: OE long 'ǣ' /æː/ shifts to ME long 'ē' /ɛː/, spelled 'e'")
    print(f"   - Spelling: OE 'hw' becomes 'wh'")
    print(f"   - Ending: The OE ending '-eþ' becomes ME '-eth'")
    print(f"   - The resulting ME infinitive is: {me_infinitive}")
    print(f"   - The final 3rd person singular present form is: {me_3sg_present}")
    print("-" * 20)

    # --- Final Derivation Summary ---
    print("Final Derivation Chain:")
    final_equation = (f"PIE: {pie_causative} -> PGmc: {pgmc_infinitive} -> "
                      f"OE: {oe_3sg_present} -> ME: {me_3sg_present}")
    print(final_equation)


if __name__ == '__main__':
    trace_word_derivation()
    final_answer = "whereth"
    print(f"\nThus, the predicted third person singular present verbal form in Middle English is '{final_answer}'.")
    # The final answer in the required format will be printed at the very end.
    print(f'<<<{final_answer}>>>')
