def solve_greek_puzzle():
    """
    This function identifies non-classical Greek words in two passages,
    provides their classical counterparts, and the historical period of the
    non-classical form. It then prints the result in the specified format.
    """
    # Passage 1 Analysis
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "KoineDemotic"

    # Passage 2 Analysis
    W2 = "ὄρνιθα"
    C2 = "ὄρνιν"
    P2 = "Homeric"

    # Combine the results into the final format
    result = f"{W1},{C1},{P1},{W2},{C2},{P2}"
    
    # Print the final result string
    print(result)

solve_greek_puzzle()