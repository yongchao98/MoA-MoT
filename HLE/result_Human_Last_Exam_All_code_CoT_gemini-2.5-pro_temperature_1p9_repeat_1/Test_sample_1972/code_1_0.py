def find_borges_reference():
    """
    Identifies the novel and author based on Jorge Luis Borges's praise.
    """
    author = "Juan Carlos Onetti"
    novel = "El astillero (The Shipyard)"
    quote_clue = "the intensity of a tiger and the variety that a chess duel can achieve"
    author_clue = "a continuator (and simplifier) of Faulkner"

    print(f"The novel Jorge Luis Borges referred to is '{novel}' by {author}.")
    print(f"\nBorges praised the novel as having \"{quote_clue}\" and described its author, {author}, as \"{author_clue}.\"")

find_borges_reference()