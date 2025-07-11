def solve_greek_puzzle():
    """
    This function identifies non-classical Greek words in two passages,
    provides their classical equivalents, and the periods of the non-classical forms.
    """
    
    # Analysis of Passage 1
    W1 = "εἶμαι"  # Non-classical word: "I am"
    C1 = "εἰμί"  # Correct classical form
    P1 = "KoineDemotic" # Period(s) of the non-classical word
    
    # Analysis of Passage 2
    W2 = "ὄρνιθα" # Non-classical word: "bird" (accusative singular)
    C2 = "ὄρνιν"  # Correct classical form
    P2 = "KoineDemotic" # Period(s) of the non-classical word
    
    # Combine into the final format
    result = f"{W1},{C1},{P1},{W2},{C2},{P2}"
    
    print(result)

solve_greek_puzzle()