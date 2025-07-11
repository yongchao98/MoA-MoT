def solve():
    """
    This function identifies non-classical words in two Greek passages,
    provides their classical equivalents and the historical periods of the
    non-classical forms, and prints them in the specified format.
    """
    # Passage 1 analysis
    w1 = "εἶμαι"
    c1 = "εἰμί"
    p1 = "KoineByzantineDemotic"

    # Passage 2 analysis
    w2 = "ὄρνιθα"
    c2 = "ὄρνιν"
    p2 = "Demotic"

    # Combine and format the results
    result = f"{w1},{c1},{p1},{w2},{c2},{p2}"
    print(result)

solve()