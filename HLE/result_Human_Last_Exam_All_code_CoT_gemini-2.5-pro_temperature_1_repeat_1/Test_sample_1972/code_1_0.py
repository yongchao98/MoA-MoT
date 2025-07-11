def find_borges_reference():
    """
    Identifies and prints the novel and author praised by Jorge Luis Borges.
    """
    # The author and novel identified from Borges's praise.
    # Borges wrote a prologue for an edition of this author's work.
    author = "Juan Rulfo"
    novel = "Pedro Páramo"
    
    # Constructing the answer sentence.
    answer = (
        f"The novel Jorge Luis Borges was referring to is '{novel}' by {author}.\n"
        f"Borges praised it in a prologue, describing it as having 'the intensity of a tiger and the variety that a chess duel can achieve,' "
        f"and called its author a continuator (and simplifier) of Faulkner."
    )
    
    print(answer)

find_borges_reference()