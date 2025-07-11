def solve_language_puzzle():
    """
    This function identifies the two languages based on the given clues
    and prints the result.
    """
    # Clues for language 'a':
    # - No 'k', 'w'
    # - Has 'à'
    # - Still in use / officially recognized
    # Analysis: Italian fits perfectly. The standard alphabet has 21 letters,
    # excluding 'k' and 'w'. The letter 'à' is common.
    language_a = "Italian"

    # Clues for language 'b':
    # - Letter combination "ggj" is widely used
    # - Letter combination "skt" is widely used
    # - Still in use / officially recognized
    # Analysis: The 'ggj' trigraph is characteristic of Faroese (e.g., 'leggja').
    # The 'skt' combination is also common in Faroese (e.g., 'føroyskt').
    language_b = "Faroese"

    print(f"Language a is: {language_a}")
    print(f"Language b is: {language_b}")

solve_language_puzzle()