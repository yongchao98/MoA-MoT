def solve_greek_puzzle():
    """
    This function identifies non-classical words in two Greek passages,
    provides their classical corrections and the periods they belong to.
    """
    # Analysis of Passage 1:
    # The word εἶμαι is the 1st person singular of 'to be'.
    # This form is characteristic of later Greek.
    W1 = "εἶμαι"
    # The correct Classical Attic form is εἰμί.
    C1 = "εἰμί"
    # The form εἶμαι developed in Koine Greek and became standard in
    # Byzantine and Demotic (Modern) Greek.
    P1 = "KoineByzantineDemotic"

    # Analysis of Passage 2:
    # The word ὄρνιθα is the accusative singular of ὄρνις ('bird').
    # This declension is a later simplification.
    W2 = "ὄρνιθα"
    # The correct Classical Attic accusative singular is ὄρνιν.
    C2 = "ὄρνιν"
    # This change from the 3rd declension to a 1st declension pattern
    # for the accusative is found in Koine and Demotic Greek.
    P2 = "KoineDemotic"

    # Format the final answer as requested.
    final_answer = f"{W1},{C1},{P1},{W2},{C2},{P2}"
    print(final_answer)

solve_greek_puzzle()