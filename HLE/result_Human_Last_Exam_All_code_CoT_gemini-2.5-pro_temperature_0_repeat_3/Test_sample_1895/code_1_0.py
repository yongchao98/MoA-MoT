def solve_semantic_transparency():
    """
    This function models the semantic transparency of languages based on the
    theory of Müller-Gotama (1994) and general linguistic principles.

    Semantic transparency refers to how easily the meaning of a complex word
    can be deduced from its parts.
    - Russian: High transparency due to a rich and regular derivational system.
    - German: High transparency due to productive compounding.
    - Old English: More transparent than Modern English, with clearer Germanic roots and compounds.
    - Modern English: Less transparent due to heavy borrowing of opaque words from Latin and French.

    We assign arbitrary scores to represent this ranking.
    """
    language_transparency_scores = {
        'Modern English': 40,
        'Old English': 70,
        'Russian': 90,
        'German': 85
    }

    # Sort the languages by their transparency score in descending order
    sorted_languages = sorted(language_transparency_scores.items(), key=lambda item: item[1], reverse=True)

    # Extract the names of the languages in the sorted order
    sorted_language_names = [lang[0] for lang in sorted_languages]

    # Print the final sequence
    print('>'.join(sorted_language_names))

solve_semantic_transparency()