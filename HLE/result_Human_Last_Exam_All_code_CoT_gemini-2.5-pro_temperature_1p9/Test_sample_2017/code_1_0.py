def solve_greek_puzzle():
    """
    This function identifies non-classical words in two Greek passages,
    provides their classical equivalents, and the period they belong to.
    """

    # Analysis of Passage 1
    # The word εἶμαι is the first-person singular "I am".
    # In Classical Greek, the form is εἰμί.
    # The form εἶμαι is found in Koine Greek and is the standard in Demotic (Modern) Greek.
    w1 = "εἶμαι"
    c1 = "εἰμί"
    p1 = "KoineDemotic"

    # Analysis of Passage 2
    # The word νηός is the genitive singular of ναῦς (ship).
    # This form is characteristic of Homeric and Ionic Greek.
    # The standard Classical Attic prose form is νεώς.
    w2 = "νηός"
    c2 = "νεώς"
    p2 = "Homeric"

    # Format the final output string
    result = f"{w1},{c1},{p1},{w2},{c2},{p2}"
    print(result)

solve_greek_puzzle()