def solve_greek_puzzle():
    """
    This function identifies anachronistic words in two Greek passages,
    provides their classical equivalents, and the periods they belong to.
    """
    
    # Analysis results for Passage 1
    # W1 is the non-classical word, C1 is the correct classical version, P1 is the period.
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "KoineDemotic"

    # Analysis results for Passage 2
    # W2 is the non-classical word, C2 is the correct classical version, P2 is the period.
    W2 = "ὄρνιθα"
    C2 = "ὄρνιν"
    P2 = "KoineDemotic"

    # Combine all parts into the final required format.
    result = f"{W1},{C1},{P1},{W2},{C2},{P2}"
    
    # Print the final result string.
    print(result)

solve_greek_puzzle()