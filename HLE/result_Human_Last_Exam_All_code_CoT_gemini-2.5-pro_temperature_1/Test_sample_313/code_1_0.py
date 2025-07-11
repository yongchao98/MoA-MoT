def solve_language_puzzle():
    """
    This script identifies two languages, 'a' and 'b', based on a set of orthographic clues.
    """

    # --- Language 'a' Identification ---
    # Clues:
    # 1. No letters 'k' or 'w' in the native orthography.
    # 2. Contains the letter 'à'.
    # 3. Is an officially recognized language used today.
    #
    # Reasoning: The combination of having 'à' while lacking native 'k' and 'w'
    # points strongly to Italian. The standard Italian alphabet has 21 letters,
    # and 'k', 'w', 'j', 'x', 'y' are only used in loanwords. It is an official
    # language of Italy and other nations.
    language_a = "Italian"

    # --- Language 'b' Identification ---
    # Clues:
    # 1. The letter combination "ggj" is very widely used.
    # 2. The letter combination "skt" is very widely used.
    # 3. Is an officially recognized language used today.
    #
    # Reasoning: The cluster "ggj" is a key identifier for Maltese, representing
    # the sound /dʒ/ in an older, but influential, orthography (now written as 'ġ').
    # The cluster "skt" is also found in Maltese, often in words borrowed from Italian.
    # Maltese is the official language of Malta and the EU.
    language_b = "Maltese"

    # --- Output the Answer ---
    print(f"Language a is: {language_a}")
    print(f"Language b is: {language_b}")

solve_language_puzzle()