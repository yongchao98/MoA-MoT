def solve_greek_puzzle():
    """
    This function identifies non-classical words in two Greek passages,
    provides their classical equivalents, and the period of the non-classical form.
    """
    # Passage 1 analysis
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "KoineDemotic"

    # Passage 2 analysis
    W2 = "ὄρνιθα"
    C2 = "ὄρνιν"
    P2 = "Koine"

    # Format and print the result
    result = f"{W1},{C1},{P1},{W2},{C2},{P2}"
    print(result)

solve_greek_puzzle()