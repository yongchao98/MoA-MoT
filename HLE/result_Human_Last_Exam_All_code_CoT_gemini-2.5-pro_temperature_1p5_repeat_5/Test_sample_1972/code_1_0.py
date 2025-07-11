def find_borges_reference():
    """
    Identifies the novel and author referenced by Jorge Luis Borges
    with the given quote and description.
    """
    novel_title_es = "La vida breve"
    novel_title_en = "A Brief Life"
    author = "Juan Carlos Onetti"
    
    quote = "the intensity of a tiger and the variety that a chess duel can achieve"
    description = "a continuator (and simplifier) of Faulkner"

    print("Jorge Luis Borges was referring to the novel:")
    print(f"- Title: {novel_title_es} ({novel_title_en})")
    print("\nHe described its author, {}, as {}.".format(author, description))

find_borges_reference()