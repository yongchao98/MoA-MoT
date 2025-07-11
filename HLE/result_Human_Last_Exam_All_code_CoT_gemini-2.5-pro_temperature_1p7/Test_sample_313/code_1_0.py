def solve_language_mystery():
    """
    This function prints the step-by-step reasoning for identifying two languages
    based on a set of clues.
    """

    print("--- Identifying Language 'a' ---")
    print("Clue 1: Still used today and is an officially recognized minority language.")
    print("Clue 2: The letters 'k' and 'w' are not in its orthography.")
    print("Clue 3: The letter 'à' is in its orthography.")
    print("\nAnalysis:")
    print("The absence of 'k' and 'w' is a strong pointer. The presence of 'à' narrows it down.")
    print("Scottish Gaelic fits perfectly: Its native alphabet does not use 'k' or 'w', it uses grave accents like 'à',")
    print("and it is an officially recognized minority language in Scotland (UK).")
    print("\nConclusion for Language a: Scottish Gaelic\n")

    print("--- Identifying Language 'b' ---")
    print("Clue 1: Still used today and is an officially recognized minority language.")
    print("Clue 2: Letter combinations 'ggj' and 'skt' are very widely used.")
    print("\nAnalysis:")
    print("The letter combination 'ggj' is very distinctive and uncommon globally.")
    print("It is a known feature of the orthography of Northern Sámi, which features consonant gradation.")
    print("Northern Sámi is also an officially recognized minority language in Norway, Sweden, and Finland.")
    print("\nConclusion for Language b: Northern Sámi\n")

    print("--- Final Answer ---")
    print("Language a is Scottish Gaelic.")
    print("Language b is Northern Sámi.")

if __name__ == '__main__':
    solve_language_mystery()