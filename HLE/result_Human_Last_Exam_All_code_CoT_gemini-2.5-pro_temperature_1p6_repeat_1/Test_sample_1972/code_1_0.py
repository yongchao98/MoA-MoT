def find_borges_reference():
    """
    This function provides the answer to the user's literary query.
    It identifies the novel and author based on Jorge Luis Borges's specific praise.
    """
    novel_title = "The Tunnel (El túnel)"
    author_name = "Ernesto Sabato"
    quote_fragment_1 = "the intensity of a tiger"
    quote_fragment_2 = "the variety that a chess duel can achieve"
    comparison_author = "Faulkner"

    print(f"The novel Jorge Luis Borges praised is '{novel_title}' by the Argentine writer {author_name}.")
    print(f"Borges wrote the prologue for the novel, in which he said it possesses \"{quote_fragment_1}\" and \"{quote_fragment_2}.\"")
    print(f"He described its author, {author_name}, as a continuator (and simplifier) of {comparison_author}, admiring how he handled dark, psychological themes in a more direct and classical manner.")

find_borges_reference()