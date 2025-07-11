def solve_greek_puzzle():
    """
    Identifies non-classical words in two Greek passages, provides their classical
    equivalents, and the periods to which the non-classical forms belong.
    """
    # Analysis of Passage 1
    w1 = "εἶμαι"
    c1 = "εἰμί"
    p1 = "KoineDemotic"

    # Analysis of Passage 2
    w2 = "ὄρνιθα"
    c2 = "ὄρνιν"
    p2 = "KoineDemotic"

    # Format the result string as per the user's request
    result = f"{w1},{c1},{p1},{w2},{c2},{p2}"
    
    print(result)

solve_greek_puzzle()