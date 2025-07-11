def find_borges_reference():
    """
    Identifies and prints the information about the novel and author
    Jorge Luis Borges referred to with the given quote.
    """
    novel_title = "Pedro Páramo"
    author_name = "Juan Rulfo"
    comparison_author = "Faulkner"
    
    quote_part_1 = "the intensity of a tiger"
    quote_part_2 = "the variety that a chess duel can achieve"
    description = f"a continuator (and simplifier) of {comparison_author}"

    print("The novel Jorge Luis Borges refers to is:")
    print(f"- Title: {novel_title}")
    print(f"- Author: {author_name}")
    print("\nBorges's description:")
    print(f"He said it has \"{quote_part_1} and {quote_part_2},\"")
    print(f"and described its author as \"{description}.\"")

find_borges_reference()