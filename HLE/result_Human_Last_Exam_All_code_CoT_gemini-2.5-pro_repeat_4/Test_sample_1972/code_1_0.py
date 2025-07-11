def find_borges_reference():
    """
    This function provides the answer to the user's literary question.
    """
    novel = "La vida breve (A Brief Life)"
    author = "Juan Carlos Onetti"
    context = (
        "Jorge Luis Borges made this reference in his prologue to the 1968 edition of the novel."
    )

    answer = (
        f"The novel Jorge Luis Borges was referring to is '{novel}' by the Uruguayan author {author}.\n\n"
        f"{context}"
    )

    print(answer)

find_borges_reference()