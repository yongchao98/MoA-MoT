def find_borges_reference():
    """
    This function stores and prints the information about the novel and author
    Jorge Luis Borges referred to with the given descriptions.
    """
    novel_title = "Pedro Páramo"
    author_name = "Juan Rulfo"
    quote_intensity = "the intensity of a tiger and the variety that a chess duel can achieve"
    quote_author = "a continuator (and simplifier) of Faulkner"

    print(f"The novel Jorge Luis Borges refers to is '{novel_title}' by the author {author_name}.")
    print("\nIn his prologue to the book, Borges describes it as having:")
    print(f"'{quote_intensity}'")
    print(f"\nHe describes the author, {author_name}, as:")
    print(f"'{quote_author}'")

if __name__ == "__main__":
    find_borges_reference()