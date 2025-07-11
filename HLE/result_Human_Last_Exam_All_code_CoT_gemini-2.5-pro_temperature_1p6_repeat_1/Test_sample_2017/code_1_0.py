def solve_greek_puzzle():
    """
    This function identifies the non-classical words in the provided Greek passages,
    their classical corrections, and the historical period of the non-classical forms.
    It then prints the result in the specified comma-separated format.
    """
    # Word 1, Classical 1, Period 1 from Passage 1
    w1 = "εἶμαι"
    c1 = "εἰμί"
    p1 = "Demotic"

    # Word 2, Classical 2, Period 2 from Passage 2
    w2 = "ὄρνιθα"
    c2 = "ὄρνιν"
    p2 = "Koine"

    # Combine all parts into the final answer string
    result = f"{w1},{c1},{p1},{w2},{c2},{p2}"
    
    print(result)

solve_greek_puzzle()