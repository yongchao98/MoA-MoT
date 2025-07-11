def solve_greek_puzzle():
    """
    This function identifies non-classical words in two Greek passages,
    provides their classical equivalents, and the periods they belong to.
    """
    # Analysis for Passage 1
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "KoineDemotic"

    # Analysis for Passage 2
    W2 = "ὄρνιθα"
    C2 = "ὄρνιν"
    P2 = "KoineDemotic"

    # Combine the results into the specified format
    result = f"{W1},{C1},{P1},{W2},{C2},{P2}"

    # Print the final result
    print(result)

solve_greek_puzzle()