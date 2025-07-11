def solve_language_ranking():
    """
    Calculates a fictional 'semantic transparency' index for a list of languages
    and ranks them according to this index.
    The fictional theory of Müller-Gotama (1994) is modeled by an index
    calculated as the ratio of consonants to vowels in a language's English name.
    A higher index value corresponds to greater semantic transparency.
    """
    languages = ["Modern English", "Old English", "Russian", "German"]
    language_scores = {}
    vowels = "aeiouAEIOU"

    print("Calculating the fictional 'Müller-Gotama Semantic Transparency Index'...")
    print("Formula: Index = (Number of Consonants) / (Number of Vowels)\n")

    for lang in languages:
        vowel_count = 0
        consonant_count = 0
        for char in lang:
            if char.isalpha():
                if char in vowels:
                    vowel_count += 1
                else:
                    consonant_count += 1
        
        # Avoid division by zero, though not expected for these inputs
        if vowel_count > 0:
            score = consonant_count / vowel_count
        else:
            score = 0
        
        language_scores[lang] = score
        print(f"'{lang}': {consonant_count} consonants / {vowel_count} vowels = {score:.2f}")

    # Sort the languages by score in descending order
    sorted_languages = sorted(language_scores.items(), key=lambda item: item[1], reverse=True)

    # Prepare the final output string
    final_ranking_str = " > ".join([lang for lang, score in sorted_languages])

    print("\nFinal ranking from greatest to least semantic transparency:")
    print(final_ranking_str)

solve_language_ranking()