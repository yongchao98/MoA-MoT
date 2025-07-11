def solve_greek_puzzle():
    """
    This function identifies non-classical Greek words, their classical versions,
    and their respective periods from two given passages.
    """
    # Passage 1 analysis
    w1 = "εἶμαι"
    c1 = "εἰμί"
    p1 = "KoineDemotic"

    # Passage 2 analysis
    w2 = "ὄρνιθα"
    c2 = "ὄρνιν"
    p2 = "Homeric"

    # Combine the results into the specified format
    result = f"{w1},{c1},{p1},{w2},{c2},{p2}"
    
    print(result)

solve_greek_puzzle()