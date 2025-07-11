def trace_word_evolution():
    """
    Traces a hypothetical PIE root through its sound changes to Middle English
    and prints the derivation.
    """
    # 1. Proto-Indo-European (PIE)
    pie_root = "*kʷeys-"
    pie_causative = "*kʷoys-éye-ti"
    print("1. The process begins with the Proto-Indo-European root and the formation of the o-grade causative:")
    print(f"   - PIE Root: {pie_root} ('to see, to heed')")
    print(f"   - PIE O-Grade Causative (3rd p. sg.): {pie_causative} ('he/she/it causes to see' -> 'he shows')")
    print("-" * 20)

    # 2. Proto-Germanic (PGmc)
    pgmc_form = "*hʷazjiþ"
    print("2. The word enters Proto-Germanic, undergoing several key changes:")
    print(f"   - PGmc Form: {pgmc_form}")
    print("   - Explanation:")
    print("     - Grimm's Law: *kʷ > *hʷ")
    print("     - Verner's Law: *s > *z (due to the accent on the PIE suffix)")
    print("     - Vowel Shift: *o > *a")
    print("     - Suffix Reduction: *-éye- > *-j- (forming a Class 1 Weak Verb)")
    print("     - Ending Change: *-ti > *-iþ")
    print("-" * 20)

    # 3. Proto-West Germanic (PWGmc)
    pwgmc_form = "*hʷarriþ"
    print("3. In Proto-West Germanic, the consonants continue to evolve:")
    print(f"   - PWGmc Form: {pwgmc_form}")
    print("   - Explanation:")
    print("     - West Germanic Gemination: *z > *zz before the following *j")
    print("     - Rhotacism: The newly doubled *zz becomes *rr")
    print("     - The *j is lost after causing the gemination")
    print("-" * 20)

    # 4. Old English (OE)
    oe_form = "hwærreþ"
    print("4. The word is inherited into Old English:")
    print(f"   - OE Form: {oe_form}")
    print("   - Explanation:")
    print("     - i-Umlaut: The root vowel *a is fronted to æ due to the 'i' in the ending *-iþ")
    print("     - Spelling & Ending: *hʷ is written 'hw', and the ending *-iþ becomes '-eþ'")
    print("-" * 20)

    # 5. Middle English (ME)
    me_form = "wherreth"
    print("5. Finally, the word evolves into its Middle English form:")
    print(f"   - ME Form: {me_form}")
    print("   - Explanation:")
    print("     - Vowel Shift: The short OE vowel 'æ' that resulted from i-umlaut becomes 'e'")
    print("     - Spelling: 'hw' is commonly respelled as 'wh'")
    print("     - Ending: The verbal ending '-eþ' is preserved as '-eth'")
    print("-" * 20)

    # Final Derivation Equation
    print("The complete derivation can be summarized as:")
    final_equation = f"{pie_causative} > {pgmc_form} > {pwgmc_form} > {oe_form} > {me_form}"
    print(final_equation)

trace_word_evolution()
<<<wherreth>>>