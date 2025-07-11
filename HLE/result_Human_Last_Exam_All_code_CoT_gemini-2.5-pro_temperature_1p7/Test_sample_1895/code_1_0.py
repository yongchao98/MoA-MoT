def solve_language_transparency():
    """
    This function demonstrates the language sequence for semantic transparency
    according to the Müller-Gotama (1994) theory, corresponding to answer B.

    The theory posits the order: German > Russian > Modern English > Old English.
    We assign scores to represent this hierarchy, where a larger value indicates
    greater semantic transparency.
    """

    # Assign scores based on the theory presented in choice B
    transparency_scores = {
        "German": 4,
        "Russian": 3,
        "Modern English": 2,
        "Old English": 1
    }

    # Sort the languages by score in descending order to build the equation
    sorted_languages = sorted(transparency_scores.items(), key=lambda item: item[1], reverse=True)

    # Build the final equation string showing the language and its score
    equation_parts = []
    for language, score in sorted_languages:
        equation_parts.append(f"{language}({score})")

    final_equation = " > ".join(equation_parts)

    print("Based on the Müller-Gotama (1994) theory, the semantic transparency order is:")
    print(final_equation)


solve_language_transparency()
<<<B>>>