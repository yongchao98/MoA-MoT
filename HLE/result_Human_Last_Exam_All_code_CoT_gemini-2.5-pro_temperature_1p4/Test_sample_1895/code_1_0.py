def solve_language_ranking():
    """
    This function solves the language ranking problem based on a fictional theory.
    It assigns scores, sorts the languages, and prints the results in a specified format.
    """
    # Step 1 & 2: Define the fictional "Müller-Gotama (1994)" Transparency Scores.
    # The theory posits that more analytic languages (less inflection) have higher semantic transparency.
    # The score is an arbitrary value from 1 to 10 representing this principle.
    transparency_scores = {
        "Modern English": 9,
        "Old English": 6,
        "Russian": 4,
        "German": 3
    }

    # Step 3: Sort the languages based on their transparency score in descending order.
    sorted_languages = sorted(transparency_scores.keys(), key=lambda lang: transparency_scores[lang], reverse=True)

    print("Based on the fictional Müller-Gotama (1994) theory, the semantic transparency scores are calculated as follows:")

    # Print the "equation" for each language as requested.
    for lang in sorted_languages:
        score = transparency_scores[lang]
        # The instruction is to "output each number in the final equation".
        # We will represent this as: Language Name Transparency Value = Score
        print(f"{lang} Transparency Value = {score}")

    # Generate the final ranked order string.
    final_order_string = " > ".join(sorted_languages)

    print("\nTherefore, the final order from greatest to least semantic transparency is:")
    print(final_order_string)

solve_language_ranking()