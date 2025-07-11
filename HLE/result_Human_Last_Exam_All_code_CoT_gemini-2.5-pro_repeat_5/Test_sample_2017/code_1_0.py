def solve_greek_puzzle():
    """
    This function identifies non-classical words in two Greek passages,
    provides their classical corrections and historical periods, and prints
    the result in the specified format.
    """
    # Analysis results
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "KoineDemotic"
    
    W2 = "ὄρνιθα"
    C2 = "ὄρνιν"
    P2 = "KoineDemotic"

    # Format the final answer string
    result = f"{W1},{C1},{P1},{W2},{C2},{P2}"
    
    print(result)

solve_greek_puzzle()