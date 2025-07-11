def solve_linguistic_derivation():
    """
    This function traces a PIE root through its historical sound changes
    to find its reflex in Middle English.
    """
    steps = [
        ("1. Proto-Indo-European (PIE) Root", "*kʷeys-", "(to see, to heed)"),
        ("2. PIE o-grade causative", "*kʷoys-éye-", "(to make see -> to show)"),
        ("3. Proto-Germanic (PGmc) form", "*hwaizijaną", "(*kʷ > *hʷ, *oy > *ai, *s > *z)"),
        ("4. Old English (OE) 3rd sg. present", "hwǣreþ", "(*ai > ǣ via i-umlaut, *z > r, + ending -eþ)"),
        ("5. Middle English (ME) 3rd sg. present", "whereth", "(hw > wh, ǣ > e /ɛː/, -eþ > -eth)")
    ]

    print("Deriving the Middle English verb for 'he shows' from PIE *kʷeys-:")
    print("-" * 60)

    final_form = ""
    # We want to print the derivation path clearly
    derivation_path = []
    for _, form, _ in steps:
        derivation_path.append(form)
    print("Full Path: " + " > ".join(derivation_path))
    print("-" * 60)


    # Print the detailed breakdown of each step
    for stage, form, notes in steps:
        print(f"{stage:<35} {form:<15} {notes}")
        if stage.startswith("5."):
            final_form = form.strip()
    
    # Final conclusion as requested
    print("-" * 60)
    print(f"The predicted third person singular present form in Middle English is '{final_form}'.")


solve_linguistic_derivation()
print("<<<whereth>>>")