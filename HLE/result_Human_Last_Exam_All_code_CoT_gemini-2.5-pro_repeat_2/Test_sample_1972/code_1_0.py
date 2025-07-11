def find_borges_reference():
    """
    This function identifies and prints the novel and author
    Borges described with the given quotes.
    """
    # Information found through literary research
    author = "Juan Carlos Onetti"
    novel = "El astillero (The Shipyard)"
    quote_1 = "the intensity of a tiger"
    quote_2 = "the variety that a chess duel can achieve"
    author_description = "a continuator (and simplifier) of Faulkner"

    # Print the findings
    print(f"Jorge Luis Borges was referring to the novel '{novel}' by the author {author}.")
    print("\nHere are the details:")
    print(f"Author: {author}")
    print(f"Novel: {novel}")
    print(f"Quote: Borges said the novel has \"{quote_1} and {quote_2}\".")
    print(f"Description: Borges described {author} as \"{author_description}\".")

find_borges_reference()