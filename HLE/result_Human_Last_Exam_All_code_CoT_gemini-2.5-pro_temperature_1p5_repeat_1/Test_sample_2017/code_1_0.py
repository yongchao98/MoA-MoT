def solve_greek_puzzle():
    """
    This function identifies non-classical words in two Greek passages,
    provides their classical equivalents, and the periods they belong to.
    """

    # Passage 1 Analysis
    # The word is 'εἶμαι'. This is the Modern Greek and Koine form.
    # The correct Classical Greek form is 'εἰμί'.
    W1 = "εἶμαι"
    C1 = "εἰμί"
    P1 = "KoineDemotic" # Common in Koine, standard in Demotic

    # Passage 2 Analysis
    # The word is 'ὄρνιθα'. This is the accusative singular based on a
    # re-declension of 'ὄρνις' to the first declension, common in later Greek.
    # The correct Classical Greek (third declension) accusative singular is 'ὄρνιν'.
    W2 = "ὄρνιθα"
    C2 = "ὄρνιν"
    P2 = "KoineDemotic" # This form emerges in Koine and is standard in Demotic

    # Combine the results into the specified format.
    result = f"{W1},{C1},{P1},{W2},{C2},{P2}"
    
    # Print the final result string.
    print(result)

solve_greek_puzzle()
<<<εἶμαι,εἰμί,KoineDemotic,ὄρνιθα,ὄρνιν,KoineDemotic>>>