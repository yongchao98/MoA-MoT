def solve_greek_puzzle():
    """
    This function provides the solution to the Greek linguistic puzzle by identifying
    non-classical words in the provided passages.
    """
    # Word 1 (W1), Correct Form 1 (C1), Period 1 (P1) from Passage 1
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "KoineByzantineDemotic"

    # Word 2 (W2), Correct Form 2 (C2), Period 2 (P2) from Passage 2
    W2 = "ὄρνιθα"
    C2 = "ὄρνιν"
    P2 = "Demotic"

    # Combine the results into the specified format
    result = f"{W1},{C1},{P1},{W2},{C2},{P2}"
    
    print(result)

solve_greek_puzzle()