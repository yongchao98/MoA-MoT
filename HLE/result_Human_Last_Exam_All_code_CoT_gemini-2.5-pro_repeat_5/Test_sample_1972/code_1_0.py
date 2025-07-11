def find_borges_reference():
    """
    Identifies and prints the novel and author Borges referred to
    based on the given quote.
    """
    novel_title = "The Big Sleep"
    author_name = "Raymond Chandler"

    print(
        f"The novel Jorge Luis Borges refers to is '{novel_title}' by {author_name}."
    )
    print(
        "\nIn a 1948 review, Borges described Chandler as a 'continuator (and simplifier) of Faulkner' "
        "and the novel's world as having 'the intensity of a tiger and the variety that a chess duel can achieve.'"
    )

if __name__ == "__main__":
    find_borges_reference()