def trace_word_development():
    """
    Traces the hypothetical development of a PIE root to a Middle English word
    and prints the derivation step by step.
    """
    # 1. Proto-Indo-European (PIE)
    pie_root = "*kʷeys-"
    pie_causative_stem = "*kʷoys-"
    print("1. Proto-Indo-European Stage")
    print(f"   - The starting root is {pie_root} ('to see, to heed').")
    print(f"   - The o-grade causative stem is {pie_causative_stem}.")
    print("-" * 20)

    # 2. Proto-Germanic (PGmc)
    pgmc_stem = "*hʷais-"
    print("2. To Proto-Germanic")
    print("   - By Grimm's Law, PIE *kʷ > PGmc *hʷ.")
    print("   - The PIE diphthong *oi becomes PGmc *ai.")
    print(f"   - Derivation: {pie_causative_stem} > {pgmc_stem}")
    print("-" * 20)

    # 3. Old English (OE)
    oe_form = "hƿǣseþ"
    print("3. To Old English (3rd Person Singular Present)")
    print("   - PGmc *ai becomes OE ā, which is then changed to ǣ by i-mutation.")
    print("   - The PGmc *hʷ becomes the digraph 'hƿ'.")
    print("   - The verb takes the 3rd person singular ending '-eþ'.")
    print(f"   - Derivation: {pgmc_stem} > {oe_form}")
    print("-" * 20)

    # 4. Middle English (ME)
    me_form = "wheseth"
    print("4. To Middle English (3rd Person Singular Present)")
    print("   - The OE 'hƿ' is respelled as 'wh'.")
    print("   - The OE vowel ǣ becomes 'e' in Middle English.")
    print("   - The ending '-eþ' becomes '-eth' in Midland dialects.")
    print(f"   - Derivation: {oe_form} > {me_form}")
    print("-" * 20)
    
    # Final Answer Summary
    print("\nFinal Derivational Path:")
    final_equation = f"{pie_causative_stem} > {pgmc_stem}- > {oe_form} > {me_form}"
    print(final_equation)
    
    # The final answer in the required format
    print(f"\nThe resulting Middle English form would be: {me_form}")
    print(f"\n<<< {me_form} >>>")

trace_word_development()