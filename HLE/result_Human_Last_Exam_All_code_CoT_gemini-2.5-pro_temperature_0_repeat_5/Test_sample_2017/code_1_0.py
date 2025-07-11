def solve_greek_puzzle():
    """
    This function identifies non-classical words in two Greek passages,
    provides their classical equivalents, and the periods they belong to.
    """
    # Passage 1 analysis results
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "KoineByzantineKatharevousaDemotic"

    # Passage 2 analysis results
    W2 = "ὄρνιθα"
    C2 = "ὄρνιν"
    P2 = "HomericKoine"

    # Combine all parts into the final answer string
    final_answer = f"{W1},{C1},{P1},{W2},{C2},{P2}"

    # Print the final answer
    print(final_answer)

solve_greek_puzzle()
print("<<<εἶμαι,εἰμί,KoineByzantineKatharevousaDemotic,ὄρνιθα,ὄρνιν,HomericKoine>>>")