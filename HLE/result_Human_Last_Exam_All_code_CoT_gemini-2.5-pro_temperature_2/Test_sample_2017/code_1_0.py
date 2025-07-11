def solve_greek_puzzle():
    """
    This function holds the solution to the Greek linguistic puzzle.
    W1: Non-classical word in Passage 1
    C1: Classical version for W1
    P1: Period for W1
    W2: Non-classical word in Passage 2
    C2: Classical version for W2
    P2: Period for W2
    """
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "KoineByzantineDemotic"
    W2 = "νηός"
    C2 = "νεώς"
    P2 = "Homeric"

    # Combine the results into the specified format.
    result = f"{W1},{C1},{P1},{W2},{C2},{P2}"
    
    print(result)

solve_greek_puzzle()