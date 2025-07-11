def solve_greek_puzzle():
    """
    This function identifies non-classical words in two Greek passages and provides the corrected versions
    and the periods the non-classical words belong to.
    """
    # Analysis results from Passage 1
    w1 = "εἶμαι"
    c1 = "εἰμί"
    p1 = "KoineDemotic"

    # Analysis results from Passage 2
    w2 = "ὄρνιθα"
    c2 = "ὄρνιν"
    p2 = "HomericKoine"

    # Combine the results into the final formatted string
    final_answer = f"{w1},{c1},{p1},{w2},{c2},{p2}"
    
    print(final_answer)

solve_greek_puzzle()