def trace_linguistic_evolution():
    """
    Traces the hypothetical evolution of a PIE root to a Middle English verb form,
    printing each step of the derivation.
    """
    pie_root = "*kʷeys-"
    target_form = "o-grade causative, 3rd person singular present"

    print(f"Task: Derive the Middle English reflex of PIE {pie_root} as an {target_form}.\n")
    print("This hypothetical word would mean 'he/she shows' (causes to see).\n")

    # Step 1: PIE to Proto-Germanic
    print("--- Step 1: Proto-Indo-European to Proto-Germanic ---")
    print("1.1: The PIE root is *kʷeys- ('to see, to heed').")
    print("1.2: A causative verb is formed with the 'o-grade' of the root vowel (*e -> *o). This gives us the stem *kʷoys-.")
    print("1.3: In Proto-Germanic, PIE causatives become Class 1 weak verbs ending in *-janą.")
    print("1.4: We apply sound changes (Grimm's Law and others):")
    print("    - PIE *kʷ -> Proto-Germanic *hw")
    print("    - PIE *o -> Proto-Germanic *a")
    print("1.5: The resulting Proto-Germanic infinitive is *hwaisijaną.")
    pgmc_verb = "*hwaisijaną"
    print(f"    Derivation: *kʷoys- + *-janą -> {pgmc_verb}\n")

    # Step 2: Proto-Germanic to Old English
    print("--- Step 2: Proto-Germanic to Old English ---")
    print(f"2.1: We start with the Proto-Germanic verb {pgmc_verb}.")
    print("2.2: The *j in the suffix causes 'i-umlaut' on the root vowel *ai.")
    print("    - Proto-Germanic *ai -> Old English *ǣ")
    print("2.3: The infinitive ending *-ijaną becomes Old English -an. The *j is lost after causing umlaut.")
    oe_infinitive = "hƿǣsan"
    print(f"2.4: The resulting Old English infinitive is '{oe_infinitive}' (using 'ƿ' for the letter wynn).")
    print("2.5: The third-person singular present ending for this class of verb is '-eþ'.")
    oe_form = "hƿǣseþ"
    print(f"    Derivation: stem hƿǣs- + ending -eþ -> {oe_form}\n")

    # Step 3: Old English to Middle English
    print("--- Step 3: Old English to Middle English ---")
    print(f"3.1: We start with the Old English form '{oe_form}'.")
    print("3.2: We apply sound and spelling changes:")
    print("    - The initial cluster hƿ /hw/ becomes spelled <wh>.")
    print("    - The long vowel *ǣ /æː/ is raised to /ɛː/, typically spelled <e>.")
    print("    - The verbal ending -eþ is preserved, spelled <-eth>.")
    me_form = "wheseth"
    print(f"3.3: The final Middle English third-person singular present form is '{me_form}'.")
    print(f"    Derivation: {oe_form} -> {me_form}\n")

    # Final Result
    print("--- Final Result ---")
    final_derivation_equation = f"PIE *kʷoys-éye-ti > PGmc. *hwaisiþi > OE. hƿǣseþ > ME. {me_form}"
    print("The complete derivation can be summarized as:")
    print(final_derivation_equation)


if __name__ == '__main__':
    trace_linguistic_evolution()
    # The final answer is the last form in the derivation.
    final_answer = "wheseth"
    print(f"\n<<<final_answer = '{final_answer}'>>>")
