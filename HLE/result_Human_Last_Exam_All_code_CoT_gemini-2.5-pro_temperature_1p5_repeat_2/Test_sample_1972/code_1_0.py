def find_borges_reference():
    """
    This function provides information about the novel and author
    referenced by Jorge Luis Borges with the given descriptions.
    """
    novel = "El astillero (The Shipyard)"
    author = "Juan Carlos Onetti"
    context = """
Jorge Luis Borges was referring to the novel "{novel}" by the Uruguayan writer {author}.

He wrote these words in his prologue for the 1961 edition of the book. In it, he praised the novel for possessing "the intensity of a tiger and the variety that a chess duel can achieve" and characterized its author, {author}, as a significant continuator (and simplifier) of the narrative techniques of William Faulkner.
    """.format(novel=novel, author=author)
    print(context)

find_borges_reference()