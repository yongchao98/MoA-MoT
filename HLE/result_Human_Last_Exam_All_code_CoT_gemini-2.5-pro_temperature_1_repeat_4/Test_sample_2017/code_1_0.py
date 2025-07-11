def solve_greek_puzzle():
    """
    This function identifies non-classical words in two Greek passages,
    provides their classical equivalents, and the periods they belong to.
    """
    # Analysis for Passage 1
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "Koine,Demotic"

    # Analysis for Passage 2
    W2 = "ὄρνιθα"
    C2 = "ὄρνιν"
    P2 = "Koine,Demotic"

    # Combine the results into the specified format
    result_string = f"{W1},{C1},{P1},{W2},{C2},{P2}"

    # Print the final result
    print(result_string)
    
    # Also provide the final answer in the required wrapper
    final_answer = f"<<<{result_string}>>>"
    print(final_answer)

solve_greek_puzzle()