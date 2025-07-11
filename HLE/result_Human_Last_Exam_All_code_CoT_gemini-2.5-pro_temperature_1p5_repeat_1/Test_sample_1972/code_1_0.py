def find_borges_reference():
    """
    This function provides the answer to the user's question about a
    novel and author referenced by Jorge Luis Borges.
    """
    novel_title = "El astillero (The Shipyard)"
    author_name = "Juan Carlos Onetti"
    quote_description = "has 'the intensity of a tiger and the variety that a chess duel can achieve'"
    author_description = "a continuator (and simplifier) of Faulkner"

    print("The novel Jorge Luis Borges refers to is:")
    print(f"- Title: {novel_title}")
    print(f"- Author: {author_name}")
    print("\nBorges described the novel as one that " + quote_description + ".")
    print("He described the author, " + author_name + ", as " + author_description + ".")

if __name__ == "__main__":
    find_borges_reference()