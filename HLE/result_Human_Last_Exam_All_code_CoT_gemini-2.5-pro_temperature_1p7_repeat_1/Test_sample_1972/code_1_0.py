def find_borges_reference():
    """
    Identifies the novel and author praised by Jorge Luis Borges
    with the given quote and description.
    """
    author = "Gabriel García Márquez"
    novel_english = "One Hundred Years of Solitude"
    novel_spanish = "Cien años de soledad"
    quote = "the intensity of a tiger and the variety that a chess duel can achieve"
    description = "a continuator (and simplifier) of Faulkner"

    print(f"The novel Jorge Luis Borges referred to is '{novel_english}' ('{novel_spanish}') by {author}.")
    print(f"Borges described it as having \"{quote},\" and its author as \"{description}.\"")

find_borges_reference()