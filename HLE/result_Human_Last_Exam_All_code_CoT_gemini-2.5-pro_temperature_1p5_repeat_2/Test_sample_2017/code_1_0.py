def solve_greek_puzzle():
    """
    This function identifies non-classical Greek words in two passages,
    provides their classical equivalents, and the period they belong to.
    """
    # Word from Passage 1, its classical form, and its period
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "KoineDemotic"

    # Word from Passage 2, its classical form, and its period
    W2 = "νηός"
    C2 = "νεώς"
    P2 = "Homeric"

    # Combine the results into the specified format
    result_string = f"{W1},{C1},{P1},{W2},{C2},{P2}"
    
    # Print the final result string
    print(result_string)

solve_greek_puzzle()
print("<<<εἶμαι,εἰμί,KoineDemotic,νηός,νεώς,Homeric>>>")