def solve_grammar_puzzle():
    """
    This function solves the Old Church Slavonic grammar puzzle provided by the user.
    """
    # The prompt asks for the form for:
    # "(3rd person singular plural, imperative, command or suggestion)"
    # We interpret the typo "singular plural" as simply "plural".
    question_description = "(3rd person plural, imperative, command or suggestion)"

    # The corresponding verb form in Old Church Slavonic is "бѫдѫтъ".
    # This is the 3rd-person plural future tense of the verb 'byti' (to be),
    # which is also used to express the imperative mood (a jussive),
    # meaning "let them be".
    answer = "бѫдѫтъ"

    # The prompt asked to "output each number in the final equation".
    # Since there are no numbers or equations, we will simply print the
    # question and its corresponding answer as a complete grammatical rule.
    print(f"{question_description} -> {answer}")

solve_grammar_puzzle()