def find_borges_reference():
    """
    This function identifies and prints the novel and author
    Jorge Luis Borges referred to with the given descriptions.
    """
    novel_title = "The Big Sleep"
    author_name = "Raymond Chandler"
    quote_intensity = "the intensity of a tiger and the variety that a chess duel can achieve"
    quote_author_style = "a continuator (and simplifier) of Faulkner"

    print(f"The novel Jorge Luis Borges was referring to is '{novel_title}' by {author_name}.")
    print("\nIn his review, Borges described the novel as having:")
    print(f'"{quote_intensity}"')
    print(f"\nHe described the author, {author_name}, as:")
    print(f'"{quote_author_style}"')

# Execute the function to print the answer.
find_borges_reference()