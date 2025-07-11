def solve_greek_puzzle():
    """
    This function identifies non-classical Greek words in two passages,
    provides their classical equivalents, and the periods of the non-classical forms.
    """
    
    # Analysis of Passage 1:
    # The word is εἶμαι (eîmai). In classical Attic Greek, the first-person singular
    # present indicative of the verb "to be" is εἰμί (eimí). The form εἶμαι
    # became standard in Koine Greek and is used in all subsequent forms of Greek.
    w1 = "εἶμαι"
    c1 = "εἰμί"
    p1 = "KoineByzantineKatharevousaDemotic"

    # Analysis of Passage 2:
    # The word is νηός (nēós), the genitive singular of ναῦς (naûs, "ship").
    # While used in poetry, the standard classical Attic prose form is νεώς (neṓs).
    # The form νηός is characteristic of Homeric and Ionic Greek.
    w2 = "νηός"
    c2 = "νεώς"
    p2 = "Homeric"

    # Format the final answer string
    result = f"{w1},{c1},{p1},{w2},{c2},{p2}"
    
    print(result)

solve_greek_puzzle()