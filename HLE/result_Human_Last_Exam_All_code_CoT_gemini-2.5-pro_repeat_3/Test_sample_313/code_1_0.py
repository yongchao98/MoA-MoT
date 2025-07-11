def solve_language_puzzle():
    """
    This function identifies the two languages based on the given clues
    and prints the result.
    """
    # Language a: No 'k' or 'w', but has 'à'.
    # This points strongly to Italian. The standard Italian alphabet has 21 letters,
    # and does not include 'k' or 'w' for native words. The grave accent 'à' is common.
    language_a = "Italian"

    # Language b: Widely uses "ggj" and "skt".
    # The "ggj" cluster (related to the Maltese letter 'ġ') and the "skt" cluster
    # are characteristic features of the Maltese language.
    language_b = "Maltese"

    print(f"Language a is: {language_a}")
    print(f"Language b is: {language_b}")

solve_language_puzzle()