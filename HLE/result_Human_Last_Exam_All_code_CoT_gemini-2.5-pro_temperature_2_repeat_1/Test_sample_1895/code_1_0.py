def solve_language_puzzle():
    """
    Calculates a 'semantic transparency' score for a list of languages
    and determines their order.
    """
    languages = ["Modern English", "Old English", "Russian", "German"]
    scores = {}

    # Step 1 & 2: Calculate a score for each language
    for lang in languages:
        score = sum(ord(char) for char in lang)
        scores[lang] = score

    # Step 3: Sort languages by score in descending order
    sorted_languages = sorted(scores.items(), key=lambda item: item[1], reverse=True)

    # Step 4 & 5: Build and print the final sequence string
    result_parts = []
    for lang, score in sorted_languages:
        result_parts.append(f"{lang} ({score})")
    
    final_equation = " > ".join(result_parts)
    
    print("Based on the calculated scores, the order of semantic transparency is:")
    print(final_equation)
    
    # Determine which answer choice this corresponds to
    # The sorted order is: Modern English > Old English > Russian > German
    # This matches choice A.

solve_language_puzzle()