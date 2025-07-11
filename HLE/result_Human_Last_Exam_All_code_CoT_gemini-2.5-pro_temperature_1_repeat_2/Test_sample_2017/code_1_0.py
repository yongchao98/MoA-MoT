def solve_greek_puzzle():
    """
    This function identifies non-classical words in two Greek passages,
    provides their classical equivalents, and the periods of the non-classical forms.
    """
    # Word from Passage 1, its classical version, and its period
    w1 = "εἶμαι"
    c1 = "εἰμί"
    p1 = "KoineDemotic"

    # Word from Passage 2, its classical version, and its period
    w2 = "ὄρνιθα"
    c2 = "ὄρνιν"
    p2 = "KoineDemotic"

    # Combine the results into the specified format
    result = f"{w1},{c1},{p1},{w2},{c2},{p2}"

    # Print the final answer in the required format
    print(f"<<<{result}>>>")

solve_greek_puzzle()