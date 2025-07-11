def find_borges_reference():
    """
    This function stores and prints the answer to the literary question
    about a quote from Jorge Luis Borges.
    """
    novel = "El astillero (The Shipyard)"
    author = "Juan Carlos Onetti"
    quote_description = "the intensity of a tiger and the variety that a chess duel can achieve"
    author_description = "a continuator (and simplifier) of Faulkner"

    print(f"The novel Jorge Luis Borges refers to is '{novel}' by the author {author}.")
    print("\nBorges described it as having:")
    print(f'"{quote_description}"')
    print(f"\nHe described the author, {author}, as:")
    print(f'"{author_description}"')

# Execute the function to print the answer.
find_borges_reference()