def solve_language_puzzle():
    """
    This function identifies two languages based on a set of linguistic clues
    and prints the result.
    """

    # Clues for Language a:
    # 1. Still in use and officially recognized (at least as a minority language).
    # 2. Orthography lacks 'k' and 'w'.
    # 3. Orthography includes 'à'.
    # Analysis: These clues strongly point to Scottish Gaelic. Its traditional
    # alphabet omits 'k' and 'w', it uses the grave accent ('à'), and it is an
    # officially recognized minority language in Scotland.
    language_a = "Scottish Gaelic"

    # Clues for Language b:
    # 1. Still in use and officially recognized.
    # 2. The letter combinations "ggj" and "skt" are widely used.
    # Analysis: The "ggj" trigraph is highly distinctive and points directly
    # to Maltese (e.g., in the word 'ħeggja'). The "skt" combination is also
    # common (e.g., in 'skutella'). Maltese is a national language of Malta.
    language_b = "Maltese"

    print(f"Language a is: {language_a}")
    print(f"Language b is: {language_b}")

solve_language_puzzle()