def solve_linguistic_puzzle():
    """
    This script traces a hypothetical Proto-Indo-European root through its
    development into Middle English, following standard linguistic rules.
    """

    # Define the stages of the word's evolution
    steps = [
        ("Proto-Indo-European Root", "*kʷeys-", "(to see, to heed)"),
        ("1. PIE O-Grade Causative Root", "*kʷoys-", "(The vowel *e becomes *o for the causative form)"),
        ("2. Proto-Germanic Root", "*hwais-", "(Via Grimm's Law: *kʷ > *hw, *o > *a)"),
        ("3. Proto-Germanic Causative Verb", "*hwaisijaną", "(Class 1 weak verb suffix *-janą is added)"),
        ("4. Old English Infinitive", "hǣsan", "(Via i-umlaut: *ai > ǣ, and loss of suffix elements)"),
        ("5. Old English 3rd Person Singular Present", "hǣseþ", "(The ending for 'he/she/it' is -eþ)"),
        ("6. Middle English 3rd Person Singular Present", "wheseth", "(OE hw > ME wh; OE ǣ > ME e; ending preserved)")
    ]

    print("Tracing the PIE root *kʷeys- to its hypothetical Middle English reflex.")
    print("-" * 70)
    
    final_derivation_chain = []

    for stage, form, explanation in steps:
        print(f"{stage:<40}: {form:<15} {explanation}")
        final_derivation_chain.append(form)

    print("-" * 70)
    print("The final derived form for 'he shows' in Middle English would be 'wheseth'.")
    print("\nThe full 'equation' or derivation chain is:")
    
    # Per the instructions, printing each "number" (interpreted as "form") in the final equation.
    print(" > ".join(final_derivation_chain))

solve_linguistic_puzzle()

# The final answer is the last form in our derivation.
final_answer = "wheseth"
print(f"\n<<<final_answer>>>")
