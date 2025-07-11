def solve_greek_passages():
    """
    Identifies non-classical words in Greek passages, provides the classical correction,
    and names the period of the non-classical word.
    """
    # Analysis for Passage 1
    w1 = "εἶμαι"
    c1 = "εἰμί"
    p1 = "KoineByzantineKatharevousaDemotic"

    # Analysis for Passage 2
    w2 = "ὄρνιθα"
    c2 = "ὄρνιν"
    p2 = "Homeric"

    # Combine into the final format
    result = f"{w1},{c1},{p1},{w2},{c2},{p2}"
    print(result)

solve_greek_passages()