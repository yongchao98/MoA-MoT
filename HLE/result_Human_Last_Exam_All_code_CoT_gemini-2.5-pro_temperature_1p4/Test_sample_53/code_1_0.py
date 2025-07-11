def trace_word_evolution():
    """
    This function traces the hypothetical evolution of a PIE root to a
    Middle English verbal form, printing each step of the derivation.
    """
    # 1. Proto-Indo-European (PIE)
    pie_root = "*kʷeys-"
    pie_causative = "*kʷoyséye-"
    print(f"Step 1: Proto-Indo-European")
    print(f"The root is {pie_root} ('to see, to heed').")
    print(f"The o-grade causative stem is {pie_causative} ('to cause to see').\n")

    # 2. Proto-Germanic (PGmc)
    pgmc_infinitive = "*hʷaizijaną"
    print(f"Step 2: Proto-Germanic")
    print(f"Sound changes applied: Grimm's Law (*kʷ > *hʷ), Verner's Law (*s > *z), and vowel shifts (*oy > *ai).")
    print(f"The resulting Class I weak verb infinitive is {pgmc_infinitive}.\n")

    # 3. Proto-West Germanic (PWGmc)
    pwgmc_infinitive = "*hʷāzijaną"
    print(f"Step 3: Proto-West Germanic")
    print(f"Sound change applied: i-umlaut causes the vowel *ai before a *j to become *ā.")
    print(f"The resulting infinitive is {pwgmc_infinitive}.\n")

    # 4. Old English (OE)
    oe_3sg_present = "hwǣrþ"
    print(f"Step 4: Old English")
    print(f"Sound changes applied: PWGmc *hʷ > OE hw, *ā > OE ǣ, and *z > r (rhotacism).")
    print(f"The infinitive becomes 'hwǣran'.")
    print(f"The 3rd person singular present form adds the ending -þ, resulting in {oe_3sg_present}.\n")

    # 5. Middle English (ME)
    me_3sg_present = "whereth"
    print(f"Step 5: Middle English")
    print(f"Sound changes applied: Spelling of hw > wh, vowel OE ǣ > ME e /ɛː/, ending -þ > -eth.")
    print(f"The resulting 3rd person singular present form is {me_3sg_present}.\n")
    
    # Final Summary
    print("--- Final Derivation Equation ---")
    final_equation_parts = [pie_causative, pgmc_infinitive, pwgmc_infinitive, oe_3sg_present, me_3sg_present]
    # The request "output each number in the final equation" is interpreted as outputting each historical form.
    print(" -> ".join(final_equation_parts))

trace_word_evolution()
<<<whereth>>>