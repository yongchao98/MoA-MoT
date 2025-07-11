def find_borges_reference():
    """
    This function identifies and prints the novel and author
    Jorge Luis Borges was referring to based on the user's query.
    """
    novel_title = "El astillero"
    english_title = "The Shipyard"
    author_name = "Juan Carlos Onetti"
    quoter_name = "Jorge Luis Borges"
    quote_fragment_1 = "the intensity of a tiger"
    quote_fragment_2 = "the variety that a chess duel can achieve"
    author_description = "a continuator (and simplifier) of Faulkner"

    print(f"{quoter_name} was referring to the novel:")
    print(f'"{novel_title}" (in English, "{english_title}")')
    print("\nWritten by the author:")
    print(f"{author_name}")
    print("\nBorges described the novel as having...")
    print(f'"{quote_fragment_1} and {quote_fragment_2},"')
    print("\n...and described its author as...")
    print(f'"{author_description}."')

find_borges_reference()